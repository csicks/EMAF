/***************************************************************************
 *
 * Authors:    Yuxuan Chen
 *
 * Copyright (C) 2021 Pattern Recognition and Bioinformatics Group, Shanghai Jiao Tong University
 *
 * Licensed under the GNU General Public License v3.0 (see LICENSE for details)
 *
 * All comments concerning this program package may be sent to e-mail address 'yxchen11@sjtu.edu.cn'
 ***************************************************************************/

#ifndef EM_IMAGE2D_H
#define EM_IMAGE2D_H

/** @file
 * this file contains data structure for two-dimensional image (real and complex)
 * TODO: may combine real/complex cases using std::complex (recommend) or _Complex/__complex__ (not recommend)
 * TODO: may replace get/set with () and remove [] (which will be supported in C++23)
 */

#include "array1D.h"
#include <random>
#include <chrono>
#include <vector>

#define MIN_D 10e-6 // threshold for judging whether two float/double value equal to each other

#define ORIENTATION_ROW 0    // operation in row orientation, i.e. x-axis
#define ORIENTATION_COLUMN 1 // operation in column orientation, i.e. y-axis

/// padding flags which is the same with opencv TODO: add more padding mode like reflection
#define PAD_CONSTANT 0  // flag for padding with constant value
#define PAD_REPLICATE 1 // flag for replicate padding, i.e. repeat value of the last pixel
#define PAD_WARP 2      // flag for warp padding, i.e. consider image as periodic signal

/** @brief
 * class for two-dimensional real image, whose size is fixed
 * for Cartesian coordinates, the positive direction for x-axis is down and the positive direction for y-axis is right
 * for polar coordinates, the positive direction for angle is counterclockwise
 */
template <typename T>
class imageReal
{
public:
    int shape[2]{};    // shape of image
    arrayReal<T> data; // data stored in one-dimensional array

public:
    imageReal()
    {
        shape[0] = 0;
        shape[1] = 0;
        data = arrayReal<T>();
    }

    explicit imageReal(const int s[2])
    {
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(static_cast<size_t>(shape[0]) * shape[1]);
    }

    imageReal(const imageReal<T> &img)
    {
        shape[0] = img.shape[0];
        shape[1] = img.shape[1];
        data = arrayReal<T>(img.data);
    }

    explicit imageReal(const imageReal<T> *img) : imageReal<T>(*img)
    {
    }

    imageReal(const arrayReal<T> &ary, const int s[2])
    {
        if (ary.length != static_cast<size_t>(s[0]) * s[1])
            throw baseException("Error: Data and shape are not matched!");
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(ary);
    }

    imageReal(const std::initializer_list<T> &initList, const int s[2])
    {
        if (initList.size() != static_cast<size_t>(s[0]) * s[1])
            throw baseException("Error: Data and shape are not matched!");
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(initList);
    }

    imageReal(const T *d, const int s[2])
    {
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(d, static_cast<size_t>(s[0]) * s[1]);
    }

    /// zoom values of image to range [0,255]
    void arrange()
    {
        data = data / data.max();
        data = (data - data.min()) / (data.max() - data.min());
        data = data * 255;
    }

    /// @return value at position (i,j), i.e. I[i,j]
    T get(int i, int j) const
    {
        return data[static_cast<size_t>(i) * shape[1] + j];
    }

    /// get value from (i,j) to (i+length,j+width), i.e. I[i:i+length,j:j+width]
    imageReal<T> get(int i, int j, int length, int width) const
    {
        if (i + length > shape[0] or j + width > shape[1])
            throw baseException("Error: Index out of range!");
        int s[2] = {length, width};
        imageReal<T> r = imageReal<T>(s);
        for (int m = 0; m < length; ++m)
            for (int n = 0; n < width; ++n)
                r.set(m, n, get(m + i, n + j));
        return r;
    }

    /// set I[i,j] to value
    void set(int i, int j, T value)
    {
        data[static_cast<size_t>(i) * shape[1] + j] = value;
    }

    /// set values which start from (i,j) to img, i.e. I[i:i+img.shape[0],j:j+img.shape[1]]=img
    void set(int i, int j, const imageReal<T> &img)
    {
        if (i + img.shape[0] > shape[0] or j + img.shape[1] > shape[1])
            throw baseException("Error: Index out of range!");
        for (int m = 0; m < img.shape[0]; ++m)
            for (int n = 0; n < img.shape[1]; ++n)
                set(m + i, n + j, img.get(m, n));
    }

    /// @return maximum value of image
    T max() const
    {
        return data.max();
    }

    /** @attention remember to free memory
     * @return index of maximum value in image, (i,j)
     */
    int *argmax() const
    {
        int *r = new int[2];
        size_t index = data.argmax();
        r[0] = index / shape[1];
        r[1] = index % shape[1];
        return r;
    }

    /// @return minimum value of image
    T min() const
    {
        return data.min();
    }

    /** @attention remember to free memory
     * @return index of maximum value in image, (i,j)
     */
    int *argmin() const
    {
        int *r = new int[2];
        size_t index = data.argmin();
        r[0] = index / shape[1];
        r[1] = index % shape[1];
        return r;
    }

    /// @return average value of image
    double mean() const
    {
        return data.mean();
    }

    /// @return sum value of image
    T sum() const
    {
        return data.sum();
    }

    /// @return variance of image
    double var() const
    {
        return data.var();
    }

    /// @return same image of this one which doesn't share memory
    imageReal<T> copy() const
    {
        imageReal<T> img = imageReal<T>(this);
        return img;
    }

    /// @return absolute value of image
    imageReal<T> abs() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.abs();
        return img;
    }

    /// @return logarithm_e value of image
    imageReal<T> log() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.log();
        return img;
    }

    /// @return sine value of image
    imageReal<T> sin() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.sin();
        return img;
    }

    /// @return cosine value of image
    imageReal<T> cos() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.cos();
        return img;
    }

    /// @return square root value of image
    imageReal<T> sqrt() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.sqrt();
        return img;
    }

    /** convert data type to another one
     * @example img.asType<double>()
     */
    template <typename T1>
    imageReal<T1> asType() const
    {
        imageReal<T1> img = imageReal<T1>(shape);
        img.data = data.template asType<T1>();
        return img;
    }

    /** @attention assignment value, i.e. img[i]=ary, doesn't work!
     * @return the ith row of image
     */
    arrayReal<T> operator[](int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[1]);
        memcpy(ary.data, data.data + static_cast<size_t>(index) * shape[1], sizeof(T) * shape[1]);
        return ary;
    }

    imageReal<T> &operator=(const imageReal<T> &img)
    {
        if (this == &img)
        {
            return *this;
        }
        shape[0] = img.shape[0];
        shape[1] = img.shape[1];
        data = img.data;
        return *this;
    }

    template <typename T1>
    imageReal<T> operator+(T1 value) const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data + value;
        return img;
    }

    template <typename T1>
    imageReal<T> operator-(T1 value) const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data - value;
        return img;
    }

    template <typename T1>
    imageReal<T> operator*(T1 value) const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data * value;
        return img;
    }

    template <typename T1>
    imageReal<T> operator/(T1 value) const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data / value;
        return img;
    }

    template <typename T1>
    imageReal<T> operator+(const imageReal<T1> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data + img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator-(const imageReal<T1> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data - img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator*(const imageReal<T1> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data * img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator/(const imageReal<T1> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data / img.data;
        return imgOut;
    }

    /// print image, with data type first
    void print() const
    {
        std::cout << typeid(T).name() << '\n';
        for (int i = 0; i < shape[0]; ++i)
        {
            for (int j = 0; j < shape[1]; ++j)
                std::cout << get(i, j) << ' ';
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /// flip image in row(x-axis) direction
    imageReal<T> flipX() const
    {
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < shape[0]; ++i)
            for (int j = 0; j < shape[1]; ++j)
            {
                img.set(i, j, get(shape[0] - 1 - i, j));
            }
        return img;
    }

    /// flip image in column(y-axis) direction
    imageReal<T> flipY() const
    {
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < shape[0]; ++i)
            for (int j = 0; j < shape[1]; ++j)
            {
                img.set(i, j, get(i, shape[1] - 1 - j));
            }
        return img;
    }

    /// flip image in both row(x-axis) and column(y-axis) direction
    imageReal<T> flipAll() const
    {
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < shape[0]; ++i)
            for (int j = 0; j < shape[1]; ++j)
            {
                img.set(i, j, get(shape[0] - 1 - i, shape[1] - 1 - j));
            }
        return img;
    }

    /// @return the ith row of image
    arrayReal<T> row(int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[1]);
        memcpy(ary.data, data.data + static_cast<size_t>(index) * shape[1], sizeof(T) * shape[1]);
        return ary;
    }

    /// @return the ith column of image
    arrayReal<T> column(int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[0]);
        for (int i = 0; i < shape[0]; ++i)
            ary[i] = get(i, index);
        return ary;
    }

    /** force the value range of image to be (minValue,maxValue)
     * value larger than maxValue will be set to maxValue
     * value smaller than minValue will be set to maxValue
     */
    imageReal<T> clip(T minValue, T maxValue) const
    {
        imageReal<T> img = imageReal<T>(this);
        img.data = img.data.clip(minValue, maxValue);
        return img;
    }

    /// replace every value in image with valueR
    imageReal<T> replace(T value, T valueR) const
    {
        imageReal<T> img = imageReal<T>(this);
        img.data = img.data.replace(value, valueR);
        return img;
    }
};

/** bi-linear interpolation of image at pixel (x,y), which doesn't consider image as periodic signal
 * @version 1
 */
template <typename T>
T interBiLinear(const imageReal<T> &img, double x, double y)
{
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    arrayReal<T> neighbour = arrayReal<T>();
    if (0 <= xl and xl < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    else
        neighbour.append(INFINITY);
    if (0 <= xl and xl < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    else
        neighbour.append(INFINITY);
    if (0 <= xr and xr < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    else
        neighbour.append(INFINITY);
    if (0 <= xr and xr < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    else
        neighbour.append(INFINITY);

    arrayReal<int> inf = neighbour.isNotInf();
    T result;

    if (inf.sum() == 0)
        result = 0;
    else if (inf.sum() == 1)
    {
        for (int i = 0; i < inf.length; ++i)
            if (inf[i] == 1)
            {
                result = neighbour[i];
                break;
            }
    }
    else
    {
        if (inf.equalList({1, 1, 0, 0}))
            result = neighbour[0] * (y - yl) + neighbour[1] * (yr - y);
        else if (inf.equalList({1, 0, 1, 0}))
            result = neighbour[0] * (x - xl) + neighbour[2] * (xr - x);
        else if (inf.equalList({0, 1, 0, 1}))
            result = neighbour[1] * (x - xl) + neighbour[3] * (xr - x);
        else if (inf.equalList({0, 0, 1, 1}))
            result = neighbour[2] * (y - yl) + neighbour[3] * (yr - y);
        else if (inf.equalList({1, 1, 1, 0}))
        {
            T temp = neighbour[0] * (y - yl) + neighbour[1] * (yr - y);
            result = temp * (x - xl) + neighbour[2] * (xr - x);
        }
        else if (inf.equalList({1, 1, 0, 1}))
        {
            T temp = neighbour[0] * (y - yl) + neighbour[1] * (yr - y);
            result = temp * (x - xl) + neighbour[3] * (xr - x);
        }
        else if (inf.equalList({1, 0, 1, 1}))
        {
            T temp = neighbour[2] * (y - yl) + neighbour[3] * (yr - y);
            result = neighbour[0] * (x - xl) + temp * (xr - x);
        }
        else if (inf.equalList({0, 1, 1, 1}))
        {
            T temp = neighbour[2] * (y - yl) + neighbour[3] * (yr - y);
            result = neighbour[1] * (x - xl) + temp * (xr - x);
        }
        else
        {
            T temp1 = neighbour[0] * (y - yl) + neighbour[1] * (yr - y);
            T temp2 = neighbour[2] * (y - yl) + neighbour[3] * (yr - y);
            result = temp1 * (x - xl) + temp2 * (xr - x);
        }
    }
    return result;
}

/// average interpolation of image at pixel (x,y), which only utilize four neighbor pixels (up,down,left,right)
template <typename T>
T interAvg(const imageReal<T> &img, double x, double y)
{
    int xl = std::floor(x);
    int xr = std::ceil(x);
    int yl = std::floor(y);
    int yr = std::ceil(y);
    arrayReal<double> neighbour = arrayReal<double>();
    if (0 <= xl and xl < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    if (0 <= xl and xl < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    if (0 <= xr and xr < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    if (0 <= xr and xr < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    if (neighbour.empty())
        return 0;
    T result = neighbour.sum() / static_cast<double>(neighbour.length);
    return result;
}

/// linear interpolation of two values
inline double linearValue(double value, double left, double right)
{
    return left + (right - left) * value;
}

/** bi-linear interpolation of image at pixel (x,y), which doesn't consider image as periodic signal
 * @version 2
 */
template <typename T>
T interBiLinear2(const imageReal<T> &img, double x, double y)
{
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    double fx = x - xl;
    double fy = y - yl;
    double neighbour[4];
    if (0 <= xl and xl < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour[0] = img.get(xl, yl);
    else
        neighbour[0] = 0;
    if (0 <= xr and xr < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour[1] = img.get(xr, yl);
    else
        neighbour[1] = 0;
    if (0 <= xl and xl < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour[2] = img.get(xl, yr);
    else
        neighbour[2] = 0;
    if (0 <= xr and xr < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour[3] = img.get(xr, yr);
    else
        neighbour[3] = 0;
    double temp0 = linearValue(fx, neighbour[0], neighbour[1]);
    double temp1 = linearValue(fx, neighbour[2], neighbour[3]);
    double r = linearValue(fy, temp0, temp1);
    return r;
}

/// move position (x,y) to be inside the image, i.e. position outside the image will be moved to the closest border
inline void adjustRange(const int s[2], double &x, double &y)
{
    while (x < 0)
        x += s[0] - 1;
    while (x > s[0] - 1)
        x -= s[0] - 1;
    while (y < 0)
        y += s[1] - 1;
    while (y > s[1] - 1)
        y -= s[1] - 1;
}

/// bi-linear interpolation of image at pixel (x,y), which considers image as periodic signal
template <typename T>
T interBiLinear2Loop(const imageReal<T> &img, double x, double y)
{
    adjustRange(img.shape, x, y);
    int xl = std::floor(x);
    xl = xl == img.shape[0] - 1 ? xl - 1 : xl;
    int xr = xl + 1;
    int yl = std::floor(y);
    yl = yl == img.shape[1] - 1 ? yl - 1 : yl;
    int yr = yl + 1;
    double fx = x - xl;
    double fy = y - yl;
    double neighbour[4];
    neighbour[0] = img.get(xl, yl);
    neighbour[1] = img.get(xr, yl);
    neighbour[2] = img.get(xl, yr);
    neighbour[3] = img.get(xr, yr);
    double temp0 = linearValue(fx, neighbour[0], neighbour[1]);
    double temp1 = linearValue(fx, neighbour[2], neighbour[3]);
    double r = linearValue(fy, temp0, temp1);
    return r;
}

/// apply warp operation to image according to affine matrix
template <typename T>
imageReal<T> warp(const imageReal<T> &img, double matrix[3][3])
{
    imageReal<T> imgNew = imageReal<T>(img.shape);
    double abs = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double matrix_inv[3][3] = {{matrix[1][1] / abs, -matrix[0][1] / abs,
                                (matrix[0][1] * matrix[1][2] -
                                 matrix[1][1] * matrix[0][2]) /
                                    abs},
                               {-matrix[1][0] / abs, matrix[0][0] / abs, (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) / abs},
                               {0, 0, 1}};
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            double x_in = matrix_inv[0][0] * i + matrix_inv[0][1] * j + matrix_inv[0][2];
            double y_in = matrix_inv[1][0] * i + matrix_inv[1][1] * j + matrix_inv[1][2];
            imgNew.set(i, j, interBiLinear2Loop(img, x_in, y_in));
        }
    return imgNew;
}

/// shift image by (x,y), positive value for down and right
template <typename T>
imageReal<T> shift(const imageReal<T> &img, double x, double y)
{
    double matrix[3][3] = {{1, 0, x},
                           {0, 1, y},
                           {0, 0, 1}};
    imageReal<T> imgNew = warp(img, matrix);
    return imgNew;
}

/// rotate image by angle, positive value for counterclockwise
template <typename T>
imageReal<T> rotate(const imageReal<T> &img, double angle)
{
    double angleR = -angle / 360 * 2 * M_PI;
    double rows = img.shape[0] / 2;
    double cols = img.shape[1] / 2;
    double matrix[3][3] = {{std::cos(angleR), std::sin(angleR),
                            rows - rows * std::cos(angleR) -
                                cols * std::sin(angleR)},
                           {-std::sin(angleR), std::cos(angleR), cols - cols * std::cos(angleR) + rows * std::sin(angleR)},
                           {0, 0, 1}};
    imageReal<T> imgNew = warp(img, matrix);
    return imgNew;
}

/// repeat array by l times to generate an image, mode meaning defined before
template <typename T>
imageReal<T> extend2Image(const arrayReal<T> &ary, int l, int mode)
{
    if (l <= 0)
        throw baseException("Error: Length of extension should positive!");
    if (mode == ORIENTATION_ROW)
    {
        int shape[2] = {l, ary.length};
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < l; ++i)
            for (int j = 0; j < ary.length; ++j)
                img.set(i, j, ary[j]);
        return img;
    }
    else if (mode == ORIENTATION_COLUMN)
    {
        int shape[2] = {ary.length, l};
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < ary.length; ++i)
            for (int j = 0; j < l; ++j)
                img.set(i, j, ary[i]);
        return img;
    }
    else
        throw baseException("Error: Orientation parameter should be 0 or 1!");
}

/// apply Hanning filter to image
template <typename T>
imageReal<T> hanning(const imageReal<T> &img)
{
    arrayReal<T> row = arrayReal<T>(img.shape[1]);
    arrayReal<T> col = arrayReal<T>(img.shape[0]);
    for (int i = 0; i < img.shape[1]; ++i)
        row[i] = 0.5 - 0.5 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    for (int i = 0; i < img.shape[0]; ++i)
        col[i] = 0.5 - 0.5 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    imageReal<T> imgRow = extend2Image(row, img.shape[0], ORIENTATION_ROW);
    imageReal<T> imgCol = extend2Image(col, img.shape[1], ORIENTATION_COLUMN);
    imageReal<T> mask = imgRow * imgCol;
    return img * mask;
}

/// apply Hamming filter to image
template <typename T>
imageReal<T> hamming(const imageReal<T> &img)
{
    arrayReal<T> row = arrayReal<T>(img.shape[1]);
    arrayReal<T> col = arrayReal<T>(img.shape[0]);
    for (int i = 0; i < img.shape[1]; ++i)
        row[i] = 0.54 - 0.46 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    for (int i = 0; i < img.shape[0]; ++i)
        col[i] = 0.54 - 0.46 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    imageReal<T> imgRow = extend2Image(row, img.shape[0], ORIENTATION_ROW);
    imageReal<T> imgCol = extend2Image(col, img.shape[1], ORIENTATION_COLUMN);
    imageReal<T> mask = imgRow * imgCol;
    return img * mask;
}

/** apply Butterworth low-pass filter to image
 * @param center: center of Butterworth low-pass filter
 * @param d0: radius of Butterworth low-pass filter
 * @param n: exponent of Butterworth low-pass filter
 */
template <typename T>
imageReal<T> butterworthLow(const imageReal<T> &img, const double center[2], double d0, int n)
{
    imageReal<T> r = imageReal<T>(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            double dn = std::sqrt((i - center[0]) * (i - center[0]) + (j - center[1]) * (j - center[1]));
            double rate = 1 + std::pow(dn / d0, n);
            r.set(i, j, img.get(i, j) / rate);
        }
    return r;
}

/** apply Butterworth high-pass filter to image
 * @param center: center of Butterworth high-pass filter
 * @param d0: radius of Butterworth high-pass filter
 * @param n: exponent of Butterworth high-pass filter
 */
template <typename T>
imageReal<T> butterworthHigh(const imageReal<T> &img, const double center[2], double d0, int n)
{
    imageReal<T> r = imageReal<T>(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            double dn = std::sqrt((i - center[0]) * (i - center[0]) + (j - center[1]) * (j - center[1]));
            if (dn == 0)
                r.set(i, j, 0);
            else
            {
                double rate = 1 + std::pow(d0 / dn, n);
                r.set(i, j, img.get(i, j) / rate);
            }
        }
    return r;
}

/** compute moment of image
 * @code
 * $\sum_{i=0}^M\sum_{j=0}^N i^p*j^q*I(i,j)$
 * @endcode
 */
template <typename T>
T moment(const imageReal<T> &img, int p, int q, int width = -1)
{
    T s = 0;
    int centerX = img.shape[0] / 2;
    int centerY = img.shape[1] / 2;
    for (int i = centerX - width; i < centerX + width; ++i)
        for (int j = centerY - width; j < centerY + width; ++j)
        {
            s += img.get(i, j) * std::pow(i, p) * std::pow(j, q);
        }
    return s;
}

/** compute center moment of image
 * @code
 * $\sum_{i=0}^M\sum_{j=0}^N (i-x_m)^p*(j-y_m)^q*I(i,j)$
 * @endcode
 */
template <typename T>
T centerMoment(const imageReal<T> &img, int p, int q, T xM, T yM, int width = -1)
{
    T s = 0;
    int centerX = img.shape[0] / 2;
    int centerY = img.shape[1] / 2;
    for (int i = centerX - width; i < centerX + width; ++i)
        for (int j = centerY - width; j < centerY + width; ++j)
        {
            s += img.get(i, j) * std::pow(i - xM, p) * std::pow(j - yM, q);
        }
    return s;
}

/// compute direction of image according to its moment
template <typename T>
double momentDirection(const imageReal<T> &img, int width = -1)
{
    if (width == -1)
        width = img.shape[0] / 2;
    T m00 = moment(img, 0, 0, width);
    T m01 = moment(img, 0, 1, width);
    T m10 = moment(img, 1, 0, width);
    T xM = m10 / m00;
    T yM = m01 / m00;
    T m11 = centerMoment(img, 1, 1, xM, yM, width);
    T m20 = centerMoment(img, 2, 0, xM, yM, width);
    T m02 = centerMoment(img, 0, 2, xM, yM, width);
    double angle = -std::atan2(2 * m11, m02 - m20) / 2 / 2 / M_PI * 360;
    return angle;
}

/// find eight nearest neighbour of position (x,y) in counterclockwise order, which starts from the upper left one
template <typename T>
arrayReal<T> findNeighbour(const imageReal<T> &img, int x, int y)
{
    int xl = x - 1;
    int xr = x + 1;
    int yl = y - 1;
    int yr = y + 1;
    arrayReal<T> neighbour = arrayReal<T>();
    if (0 <= xl and xl < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    if (0 <= xl and xl < img.shape[0])
        neighbour.append(img.get(xl, y));
    if (0 <= xl and xl < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    if (0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(x, yr));
    if (0 <= xr and xr < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    if (0 <= xr and xr < img.shape[0])
        neighbour.append(img.get(xr, y));
    if (0 <= xr and xr < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    if (0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(x, yl));
    return neighbour;
}

/// find connected domain in image
template <typename T>
imageReal<int> connectDomainTwoPass(const imageReal<T> &img)
{
    int label = 1;
    std::vector<arrayReal<int>> equalLabels;
    imageReal<int> r = imageReal<int>(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
            if (img.get(i, j) == 1)
            {
                arrayReal<int> neighbour = findNeighbour(r, i, j);
                int m1 = neighbour.max();
                if (m1 == 0)
                {
                    r.set(i, j, label);
                    ++label;
                }
                else
                {
                    neighbour = neighbour.replace(0, neighbour.max());
                    r.set(i, j, neighbour.min());
                    if (neighbour.min() != neighbour.max())
                    {
                        if (equalLabels.empty())
                        {
                            equalLabels.push_back(neighbour.unique());
                        }
                        else
                        {
                            bool flag = true;
                            for (auto &equalLabel : equalLabels)
                                if (equalLabel.intersectionA(neighbour).length != 0)
                                {
                                    equalLabel = equalLabel.unionA(neighbour);
                                    flag = false;
                                    break;
                                }
                            if (flag)
                                equalLabels.push_back(neighbour.unique());
                        }
                    }
                }
            }

    for (int i = 0; i < r.shape[0]; ++i)
        for (int j = 0; j < r.shape[1]; ++j)
            if (r.get(i, j) != 0)
                for (auto &equalLabel : equalLabels)
                    if (equalLabel.exist(r.get(i, j)) != 0)
                    {
                        r.set(i, j, equalLabel.min());
                        break;
                    }

    return r;
}

/// find largest connected domain in image
template <typename T>
imageReal<int> largestDomain(const imageReal<T> &img)
{
    imageReal<int> domain = connectDomainTwoPass(img);
    arrayReal<int> flatten = domain.data.sort();
    arrayReal<int> valueGroup;
    arrayReal<int> countGroup;
    int c = -1;
    for (int i = 0; i < flatten.length; ++i)
    {
        if (flatten[i] != 0)
        {
            if (valueGroup.exist(flatten[i]) == 0)
            {
                valueGroup.append(flatten[i]);
                countGroup.append(1);
                ++c;
            }
            else
            {
                ++countGroup[c];
            }
        }
    }
    size_t index = countGroup.argmax();
    int valueMost = valueGroup[index];
    for (int i = 0; i < domain.shape[0]; ++i)
        for (int j = 0; j < domain.shape[1]; ++j)
        {
            if (domain.get(i, j) == valueMost)
                domain.set(i, j, 1);
            else
                domain.set(i, j, 0);
        }
    return domain;
}

/** compute correntropy between two images
 * @code
 * $\sum_{i=0}^M\sum_{j=0}^N\exp(-(I_1(i,j)-I_2(i,j))^2)$
 * @endcode
 */
template <typename T>
T correntropy(const imageReal<T> &img1, const imageReal<T> &img2)
{
    if (img1.shape[0] != img2.shape[0] or img1.shape[1] != img2.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    T value = 0;
    for (int i = 0; i < img1.shape[0]; ++i)
        for (int j = 0; j < img1.shape[1]; ++j)
            value += std::exp(-std::pow(img1.get(i, j) - img2.get(i, j), 2));
    return value / img1.shape[0] / img1.shape[1];
}

/** convert image to polar coordinates
 * @implements map each pixel in polar coordinates to Cartesian coordinates
 */
template <typename T>
imageReal<T> polarBack(const imageReal<T> &img, bool logFlag = false)
{
    double centerX = img.shape[0] / 2;
    double centerY = img.shape[1] / 2;
    imageReal<T> imgPolar = imageReal<T>(img.shape);
    double maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            double angle = static_cast<double>(i) / img.shape[0] * 360;
            double radius;
            if (logFlag)
                radius = std::exp(static_cast<double>(j) / img.shape[1] * std::log(maxRadius));
            else
                radius = static_cast<double>(j) / img.shape[1] * maxRadius;
            double x = -radius * std::sin(angle / 360 * M_PI * 2) + centerX;
            double y = radius * std::cos(angle / 360 * M_PI * 2) + centerY;
            imgPolar.set(i, j, interBiLinear2Loop(img, x, y));
        }
    return imgPolar;
}

/// apply interpolation for all pixels equaling to input value
template <typename T>
imageReal<T> valueInter(const imageReal<T> &img, T value)
{
    imageReal<T> imgNew = img.copy();
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
            if (std::abs(img.get(i, j) - value) < MIN_D)
            {
                T m = interBiLinear2Loop(img, i, j);
                imgNew.set(i, j, m);
            }
    return imgNew;
}

/** convert image to polar coordinates
 * @implements map each pixel in polar Cartesian coordinates to polar coordinates, and then apply interpolation
 * @deprecated
 */
template <typename T>
imageReal<T> polarFront(const imageReal<T> &img, bool logFlag = false)
{
    double centerX = img.shape[0] / 2;
    double centerY = img.shape[1] / 2;
    imageReal<T> imgPolar = imageReal<T>(img.shape);
    double maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            double radius = std::sqrt(std::pow(i - centerX, 2) + std::pow(j - centerY, 2));
            double yD = centerX - i;
            double xD = j - centerY;
            double angle = std::atan2(yD, xD) / 2 / M_PI * 360;
            angle = angle < 0 ? angle + 360 : angle;
            int radiusIndex;
            if (logFlag)
                radiusIndex = radius == 0 ? 0 : std::round(std::log(radius) / std::log(maxRadius) * (img.shape[1] - 1));
            else
                radiusIndex = std::round(radius / maxRadius * (img.shape[1] - 1));
            int angIndex = std::round(angle / 360 * (img.shape[0] - 1));
            imgPolar.set(angIndex, radiusIndex, img.get(i, j));
        }
    imgPolar = valueInter(imgPolar, static_cast<T>(0));
    return imgPolar;
}

/** pad image using different mode, with mode defined before
 * @param top: add rows at the top of image
 * @param bottom: add rows at the bottom of image
 * @param left: add columns at the left of image
 * @param right: add rows at right top of image
 * @param value: constant value if mode PAD_CONSTANT is used
 */
template <typename T>
imageReal<T> padding(const imageReal<T> &img, int top, int bottom, int left, int right, int mode, int value = 0)
{
    int s[2] = {img.shape[0] + top + bottom, img.shape[1] + left + right};
    imageReal<T> r = imageReal<T>(s);
    switch (mode)
    {
    case PAD_CONSTANT:
        for (int j = 0; j < r.shape[1]; ++j)
        {
            for (int i = 0; i < top; ++i)
                r.set(i, j, value);
            for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
                r.set(i, j, value);
        }
        for (int i = top; i < r.shape[0] - bottom; ++i)
        {
            for (int j = 0; j < left; ++j)
                r.set(i, j, value);
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, value);
        }
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = left; j < r.shape[1] - right; ++j)
                r.set(i, j, img.get(i - top, j - left));
        break;
    case PAD_REPLICATE:
        for (int i = 0; i < top; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(0, 0));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(i - top, 0));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(i - top, img.shape[1] - 1));
        for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(img.shape[0] - 1, 0));
        for (int j = left; j < r.shape[1] - right; ++j)
            for (int i = 0; i < top; ++i)
                r.set(i, j, img.get(0, j - left));
        for (int j = left; j < r.shape[1] - right; ++j)
            for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
                r.set(i, j, img.get(img.shape[0] - 1, j - left));
        for (int i = 0; i < top; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(0, img.shape[1] - 1));
        for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(img.shape[0] - 1, img.shape[1] - 1));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = left; j < r.shape[1] - right; ++j)
                r.set(i, j, img.get(i - top, j - left));
        break;
    case PAD_WARP:
        for (int i = 0; i < top; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(img.shape[0] - top + i, img.shape[1] - left + j));
        for (int i = 0; i < top; ++i)
            for (int j = left; j < r.shape[1] - right; ++j)
                r.set(i, j, img.get(img.shape[0] - top + i, j - left));
        for (int i = 0; i < top; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(img.shape[0] - top + i, j - left - img.shape[1]));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(i - top, img.shape[1] - left + j));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(i - top, j - left - img.shape[1]));
        for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
            for (int j = 0; j < left; ++j)
                r.set(i, j, img.get(i - top - img.shape[0], img.shape[1] - left + j));
        for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
            for (int j = left; j < r.shape[1] - right; ++j)
                r.set(i, j, img.get(i - top - img.shape[0], j - left));
        for (int i = r.shape[0] - bottom; i < r.shape[0]; ++i)
            for (int j = r.shape[1] - right; j < r.shape[1]; ++j)
                r.set(i, j, img.get(i - top - img.shape[0], j - left - img.shape[1]));
        for (int i = top; i < r.shape[0] - bottom; ++i)
            for (int j = left; j < r.shape[1] - right; ++j)
                r.set(i, j, img.get(i - top, j - left));
        break;
    default:
        throw baseException("Error: Input padding mode is not supported!");
    }
    return r;
}

/** @brief
 * class for two-dimensional complex image
 * same direction as real image
 * functions/variables appeared before in imageReal are ignored
 */
class imageComplex
{
public:
    int shape[2]{};
    arrayComplex data;

public:
    imageComplex();

    explicit imageComplex(const int s[2]);

    imageComplex(const imageComplex &img);

    explicit imageComplex(const imageComplex *img);

    imageComplex(const arrayComplex &d, const int s[2]);

    imageComplex(const double d[][2], const int s[2]);

    /// @attention remember to free memory
    double *get(int i, int j) const;

    imageComplex get(int i, int j, int length, int width) const;

    void set(int i, int j, const double value[2]);

    void set(int i, int j, const imageComplex &img);

    /// @attention remember to free memory
    double *mean() const;

    imageComplex copy() const;

    imageReal<double> abs() const;

    /// @return real part of complex image
    imageReal<double> real() const;

    /// @return imaginary part of complex image
    imageReal<double> imag() const;

    /// @return complex angle of complex image
    imageReal<double> angle() const;

    arrayComplex operator[](int index) const;

    imageComplex &operator=(const imageComplex &img);

    imageComplex operator+(double value) const;

    imageComplex operator-(double value) const;

    imageComplex operator*(double value) const;

    imageComplex operator/(double value) const;

    imageComplex operator+(const double value[2]) const;

    imageComplex operator-(const double value[2]) const;

    imageComplex operator*(const double value[2]) const;

    imageComplex operator/(const double value[2]) const;

    void print() const;

    /// @return conjugate of complex image
    imageComplex conj() const;

    imageComplex operator+(const imageComplex &img) const;

    imageComplex operator-(const imageComplex &img) const;

    imageComplex operator*(const imageComplex &img) const;

    imageComplex operator/(const imageComplex &img) const;

    template <typename T>
    imageComplex operator+(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data + img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator-(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data - img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator*(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data * img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator/(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data / img.data;
        return imgOut;
    }

    imageComplex flipX() const;

    imageComplex flipY() const;

    imageComplex flipAll() const;

    arrayComplex row(int index) const;

    arrayComplex column(int index) const;
};

/// @attention remember to free memory
double *interBiLinear(const imageComplex &img, double x, double y);

/// @attention remember to free memory
double *interAvg(const imageComplex &img, double x, double y);

imageComplex warp(const imageComplex &img, double matrix[3][3]);

imageComplex shift(const imageComplex &img, double x, double y);

imageComplex rotate(const imageComplex &img, double angle);

imageComplex hanning(const imageComplex &img);

imageComplex hamming(const imageComplex &img);

imageComplex butterworthLow(const imageComplex &img, const int center[2], double d0, int n);

imageComplex butterworthHigh(const imageComplex &img, const int center[2], double d0, int n);

imageComplex polarBack(const imageComplex &img, bool logFlag = false);

arrayComplex findNeighbour(const imageComplex &img, int x, int y);

imageComplex valueInter(const imageComplex &img, const double value[2]);

imageComplex polarFront(const imageComplex &img, bool logFlag = false);

/// generate 0 or 1 round mask of radius
imageReal<int> circularMask(const int shape[2], double radius);

/// recover complex image according to its absolute value and complex angle
imageComplex bind(const imageReal<double> &abs, const imageReal<double> &angle);

#endif // EM_IMAGE2D_H
