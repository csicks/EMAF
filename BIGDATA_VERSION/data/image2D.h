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

#ifndef ALIGNMENT_IMAGE2D_H
#define ALIGNMENT_IMAGE2D_H

#include "array1D.h"
#include <random>
#include <chrono>
#include <vector>

#define MIN_D 10e-6

#define ORIENTATION_ROW 0
#define ORIENTATION_COLUMN 1

#define PAD_CONSTANT 0
#define PAD_REPLICATE 1
#define PAD_WARP 2

template <typename T>
class imageReal
{
public:
    long long shape[2]{};
    arrayReal<T> data;

public:
    imageReal()
    {
        shape[0] = 0;
        shape[1] = 0;
        data = arrayReal<T>();
    }

    explicit imageReal(const long long s[2])
    {
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(shape[0] * shape[1]);
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

    imageReal(const arrayReal<T> &d, const long long *s)
    {
        if (d.length != s[0] * s[1])
            throw baseException("Error: Data and shape are not matched!");
        shape[0] = s[0];
        shape[1] = s[1];
        data = arrayReal<T>(d);
    }

    imageReal(const std::initializer_list<T> &initList, const long long s[2])
    {
        if (initList.size() != s[0] * s[1])
            throw baseException("Error: Data and shape are not matched!");
        data = arrayReal<T>(initList);
        shape[0] = s[0];
        shape[1] = s[1];
    }

    imageReal(const T *d, const long long s[2])
    {
        data = arrayReal<T>(d, s[0] * s[1]);
        shape[0] = s[0];
        shape[1] = s[1];
    }

    void arrange()
    {
        data = data / data.max();
        data = (data - data.min()) / (data.max() - data.min());
        data = data * 255;
    }

    T get(int i, int j) const
    {
        return data[i * shape[1] + j];
    }

    imageReal<T> get(int i, int j, int height, int width) const
    {
        if (i + height > shape[0] || j + width > shape[1])
            throw baseException("Error: Index out of range!");
        long long s[2] = {height, width};
        imageReal<T> r = imageReal<T>(s);
        for (int m = 0; m < height; ++m)
            for (int n = 0; n < width; ++n)
                r.set(m, n, get(m + i, n + j));
        return r;
    }

    void set(int i, int j, T value)
    {
        data[i * shape[1] + j] = value;
    }

    void set(int i, int j, const imageReal<T> &img)
    {
        if (i + img.shape[0] > shape[0] || j + img.shape[1] > shape[1])
            throw baseException("Error: Index out of range!");
        for (int m = 0; m < img.shape[0]; ++m)
            for (int n = 0; n < img.shape[1]; ++n)
                set(m + i, n + j, img.get(m, n));
    }

    T max() const
    {
        return data.max();
    }

    int *argmax() const
    {
        int *r = new int[2];
        int index = data.argmax();
        r[0] = index / shape[1];
        r[1] = index % shape[1];
        return r;
    }

    T min() const
    {
        return data.min();
    }

    int *argmin() const
    {
        int *r = new int[2];
        int index = data.argmin();
        r[0] = index / shape[1];
        r[1] = index % shape[1];
        return r;
    }

    float mean() const
    {
        return data.mean();
    }

    T sum() const
    {
        return data.sum();
    }

    float var() const
    {
        return data.var();
    }

    imageReal<T> copy() const
    {
        imageReal<T> img = imageReal<T>(this);
        return img;
    }

    imageReal<T> abs() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.abs();
        return img;
    }

    imageReal<T> log() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.log();
        return img;
    }

    imageReal<T> sin() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.sin();
        return img;
    }

    imageReal<T> cos() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.cos();
        return img;
    }

    imageReal<T> sqrt() const
    {
        imageReal<T> img = imageReal<T>(shape);
        img.data = data.sqrt();
        return img;
    }

    template <typename T_IN>
    imageReal<T_IN> asType() const
    {
        imageReal<T_IN> img = imageReal<T_IN>(shape);
        img.data = data.template asType<T_IN>();
        return img;
    }

    arrayReal<T> operator[](int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[1]);
        memcpy(ary.data, data.data + index * shape[1], sizeof(T) * shape[1]);
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
    imageReal<T> operator+(const imageReal<T1> img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data + img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator-(const imageReal<T1> img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data - img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator*(const imageReal<T1> img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data * img.data;
        return imgOut;
    }

    template <typename T1>
    imageReal<T> operator/(const imageReal<T1> img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageReal<T> imgOut = imageReal<T>(shape);
        imgOut.data = data / img.data;
        return imgOut;
    }

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

    arrayReal<T> row(int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[1]);
        memcpy(ary.data, data.data + index * shape[1], sizeof(T) * shape[1]);
        return ary;
    }

    arrayReal<T> column(int index) const
    {
        arrayReal<T> ary = arrayReal<T>(shape[0]);
        for (int i = 0; i < shape[0]; ++i)
            ary[i] = get(i, index);
        return ary;
    }

    imageReal<T> clip(T minValue, T maxValue) const
    {
        imageReal<T> img = imageReal<T>(this);
        img.data = img.data.clip(minValue, maxValue);
        return img;
    }

    imageReal<T> replace(T value, T valueR) const
    {
        imageReal<T> img = imageReal<T>(this);
        img.data = img.data.replace(value, valueR);
        return img;
    }
};

template <typename T>
T interBiLinear(const imageReal<T> &img, float x, float y)
{
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    arrayReal<T> neighbour = arrayReal<T>();
    if (0 <= xl && xl < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    else
        neighbour.append(INFINITY);
    if (0 <= xl && xl < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    else
        neighbour.append(INFINITY);
    if (0 <= xr && xr < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    else
        neighbour.append(INFINITY);
    if (0 <= xr && xr < img.shape[0] && 0 <= yr && yr < img.shape[1])
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

template <typename T>
T interAvg(const imageReal<T> &img, float x, float y)
{
    int xl = std::floor(x);
    int xr = std::ceil(x);
    int yl = std::floor(y);
    int yr = std::ceil(y);
    arrayReal<float> neighbour = arrayReal<float>();
    if (0 <= xl && xl < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    if (0 <= xl && xl < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    if (0 <= xr && xr < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    if (0 <= xr && xr < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    if (neighbour.empty())
        return 0;
    T result = neighbour.sum() / neighbour.length;
    return result;
}

inline float linearValue(float value, float left, float right)
{
    return left + (right - left) * value;
}

template <typename T>
T interBiLinear2(const imageReal<T> &img, float x, float y)
{
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    float fx = x - xl;
    float fy = y - yl;
    float neighbour[4];
    if (0 <= xl && xl < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour[0] = img.get(xl, yl);
    else
        neighbour[0] = 0;
    if (0 <= xr && xr < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour[1] = img.get(xr, yl);
    else
        neighbour[1] = 0;
    if (0 <= xl && xl < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour[2] = img.get(xl, yr);
    else
        neighbour[2] = 0;
    if (0 <= xr && xr < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour[3] = img.get(xr, yr);
    else
        neighbour[3] = 0;
    double temp0 = linearValue(fx, neighbour[0], neighbour[1]);
    double temp1 = linearValue(fx, neighbour[2], neighbour[3]);
    double r = linearValue(fy, temp0, temp1);
    return r;
}

inline void adjustRange(const long long shape[2], float &x, float &y)
{
    while (x < 0)
        x += shape[0] - 1;
    while (x > shape[0] - 1)
        x -= shape[0] - 1;
    while (y < 0)
        y += shape[1] - 1;
    while (y > shape[1] - 1)
        y -= shape[1] - 1;
}

template <typename T>
T interBiLinear2Loop(const imageReal<T> &img, float x, float y)
{
    adjustRange(img.shape, x, y);
    int xl = std::floor(x);
    xl = xl == img.shape[0] - 1 ? xl - 1 : xl;
    int xr = xl + 1;
    int yl = std::floor(y);
    yl = yl == img.shape[1] - 1 ? yl - 1 : yl;
    int yr = yl + 1;
    float fx = x - xl;
    float fy = y - yl;
    float neighbour[4];
    neighbour[0] = img.get(xl, yl);
    neighbour[1] = img.get(xr, yl);
    neighbour[2] = img.get(xl, yr);
    neighbour[3] = img.get(xr, yr);
    double temp0 = linearValue(fx, neighbour[0], neighbour[1]);
    double temp1 = linearValue(fx, neighbour[2], neighbour[3]);
    double r = linearValue(fy, temp0, temp1);
    return r;
}

template <typename T>
imageReal<T> warp(const imageReal<T> &img, float matrix[3][3])
{
    imageReal<T> imgNew = imageReal<T>(img.shape);
    float abs = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    float matrix_inv[3][3] = {{matrix[1][1] / abs, -matrix[0][1] / abs, (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) / abs},
                              {-matrix[1][0] / abs, matrix[0][0] / abs, (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) / abs},
                              {0, 0, 1}};
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            float x_in = matrix_inv[0][0] * i + matrix_inv[0][1] * j + matrix_inv[0][2];
            float y_in = matrix_inv[1][0] * i + matrix_inv[1][1] * j + matrix_inv[1][2];
            imgNew.set(i, j, interBiLinear2Loop(img, x_in, y_in));
        }
    return imgNew;
}

template <typename T>
imageReal<T> shift(const imageReal<T> &img, float x, float y)
{
    float matrix[3][3] = {{1, 0, x},
                          {0, 1, y},
                          {0, 0, 1}};
    imageReal<T> imgNew = warp(img, matrix);
    return imgNew;
}

template <typename T>
imageReal<T> rotate(const imageReal<T> &img, float ang)
{
    float ang_r = -ang / 360 * 2 * M_PI;
    float rows = img.shape[0] / 2;
    float cols = img.shape[1] / 2;
    float matrix[3][3] = {{std::cos(ang_r), std::sin(ang_r), rows - rows * std::cos(ang_r) - cols * std::sin(ang_r)},
                          {-std::sin(ang_r), std::cos(ang_r), cols - cols * std::cos(ang_r) + rows * std::sin(ang_r)},
                          {0, 0, 1}};
    imageReal<T> imgNew = warp(img, matrix);
    return imgNew;
}

template <typename T>
imageReal<T> extend2Image(const arrayReal<T> &ary, int l, int mode)
{
    if (l <= 0)
        throw baseException("Error: Length of extension should positive!");
    if (mode == ORIENTATION_ROW)
    {
        long long shape[2] = {l, ary.length};
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < l; ++i)
            for (int j = 0; j < ary.length; ++j)
                img.set(i, j, ary[j]);
        return img;
    }
    else if (mode == ORIENTATION_COLUMN)
    {
        long long shape[2] = {ary.length, l};
        imageReal<T> img = imageReal<T>(shape);
        for (int i = 0; i < ary.length; ++i)
            for (int j = 0; j < l; ++j)
                img.set(i, j, ary[i]);
        return img;
    }
    else
        throw baseException("Error: Orientation parameter should be 0 or 1!");
}

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

template <typename T>
imageReal<T> butterworthLow(const imageReal<T> &img, const float center[2], float d0, int n)
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

template <typename T>
imageReal<T> butterworthHigh(const imageReal<T> &img, const float center[2], float d0, int n)
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

template <typename T>
float momentDirection(const imageReal<T> &img, int width = -1)
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
    float ang = -std::atan2(2 * m11, m02 - m20) / 2 / 2 / M_PI * 360;
    return ang;
}

/// eight nearest neighbour in counterclockwise order
template <typename T>
arrayReal<T> findNeighbour(const imageReal<T> &img, int x, int y)
{
    int xl = x - 1;
    int xr = x + 1;
    int yl = y - 1;
    int yr = y + 1;
    arrayReal<T> neighbour = arrayReal<T>();
    if (0 <= xl && xl < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    if (0 <= xl && xl < img.shape[0])
        neighbour.append(img.get(xl, y));
    if (0 <= xl && xl < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    if (0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(x, yr));
    if (0 <= xr && xr < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    if (0 <= xr && xr < img.shape[0])
        neighbour.append(img.get(xr, y));
    if (0 <= xr && xr < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    if (0 <= yl && yl < img.shape[1])
        neighbour.append(img.get(x, yl));
    return neighbour;
}

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
    int index = countGroup.argmax();
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

template <typename T>
T correntropy(const imageReal<T> &img1, const imageReal<T> &img2)
{
    if (img1.shape[0] != img2.shape[0] || img1.shape[1] != img2.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    T value = 0;
    for (int i = 0; i < img1.shape[0]; ++i)
        for (int j = 0; j < img1.shape[1]; ++j)
            value += std::exp(-std::pow(img1.get(i, j) - img2.get(i, j), 2));
    return value / img1.shape[0] / img1.shape[1];
}

template <typename T>
imageReal<T> polarBack(const imageReal<T> &img, bool logFlag = false)
{
    float centerX = img.shape[0] / 2;
    float centerY = img.shape[1] / 2;
    imageReal<T> imgPolar = imageReal<T>(img.shape);
    float maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            float ang = static_cast<float>(i) / img.shape[0] * 360;
            float radius;
            if (logFlag)
                radius = std::exp(static_cast<float>(j) / img.shape[1] * std::log(maxRadius));
            else
                radius = static_cast<float>(j) / img.shape[1] * maxRadius;
            float x = -radius * std::sin(ang / 360 * M_PI * 2) + centerX;
            float y = radius * std::cos(ang / 360 * M_PI * 2) + centerY;
            imgPolar.set(i, j, interBiLinear2Loop(img, x, y));
        }
    return imgPolar;
}

template <typename T>
imageReal<T> valueInter(const imageReal<T> &img, T value)
{
    imageReal<T> imgNew = img.copy();
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
            if (std::abs(img.get(i, j) - value) < MIN_D)
            {
                //                arrayReal<T> neighbour = findNeighbour(img, i, j);
                //                T m = static_cast<T>(neighbour.mean());
                T m = interBiLinear2Loop(img, i, j);
                imgNew.set(i, j, m);
            }
    return imgNew;
}

template <typename T>
imageReal<T> polarFront(const imageReal<T> &img, bool logFlag = false)
{
    float centerX = img.shape[0] / 2;
    float centerY = img.shape[1] / 2;
    imageReal<T> imgPolar = imageReal<T>(img.shape);
    float maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            float radius = std::sqrt(std::pow(i - centerX, 2) + std::pow(j - centerY, 2));
            float yD = centerX - i;
            float xD = j - centerY;
            float ang = std::atan2(yD, xD) / 2 / M_PI * 360;
            ang = ang < 0 ? ang + 360 : ang;
            int radiusIndex;
            if (logFlag)
                radiusIndex = radius == 0 ? 0 : std::round(std::log(radius) / std::log(maxRadius) * (img.shape[1] - 1));
            else
                radiusIndex = std::round(radius / maxRadius * (img.shape[1] - 1));
            int angIndex = std::round(ang / 360 * (img.shape[0] - 1));
            imgPolar.set(angIndex, radiusIndex, img.get(i, j));
        }
    imgPolar = valueInter(imgPolar, static_cast<T>(0));
    return imgPolar;
}

template <typename T>
imageReal<T> padding(const imageReal<T> &img, int top, int bottom, int left, int right, int mode, int value = 0)
{
    long long s[2] = {img.shape[0] + top + bottom, img.shape[1] + left + right};
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

class imageComplex
{
public:
    long long shape[2]{};
    arrayComplex data;

public:
    imageComplex();

    explicit imageComplex(const long long s[2]);

    imageComplex(const imageComplex &img);

    explicit imageComplex(const imageComplex *img);

    imageComplex(const arrayComplex &d, const long long *s);

    imageComplex(const double d[][2], const long long s[2]);

    double *get(int i, int j) const;

    imageComplex get(int i, int j, int height, int width) const;

    void set(int i, int j, const double value[2]);

    void set(int i, int j, const imageComplex &img);

    double *mean() const;

    imageComplex copy() const;

    imageReal<double> abs() const;

    imageReal<double> real() const;

    imageReal<double> imag() const;

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

    imageComplex conj() const;

    imageComplex operator+(const imageComplex &img) const;

    imageComplex operator-(const imageComplex &img) const;

    imageComplex operator*(const imageComplex &img) const;

    imageComplex operator/(const imageComplex &img) const;

    template <typename T>
    imageComplex operator+(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data + img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator-(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data - img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator*(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
            throw baseException("Error: Two images are not of the same size!");
        imageComplex imgOut = imageComplex(shape);
        imgOut.data = data * img.data;
        return imgOut;
    }

    template <typename T>
    imageComplex operator/(const imageReal<T> &img) const
    {
        if (shape[0] != img.shape[0] || shape[1] != img.shape[1])
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

double *interBiLinear(const imageComplex &img, float x, float y);

double *interAvg(const imageComplex &img, float x, float y);

imageComplex warp(const imageComplex &img, float matrix[3][3]);

imageComplex shift(const imageComplex &img, float x, float y);

imageComplex rotate(const imageComplex &img, float ang);

imageComplex hanning(const imageComplex &img);

imageComplex hamming(const imageComplex &img);

imageComplex butterworthLow(const imageComplex &img, const int center[2], float d0, int n);

imageComplex butterworthHigh(const imageComplex &img, const int center[2], float d0, int n);

imageComplex polarBack(const imageComplex &img, bool logFlag = false);

/// eight neighbour in counterclockwise order
arrayComplex findNeighbour(const imageComplex &img, int x, int y);

imageComplex valueInter(const imageComplex &img, const double value[2]);

imageComplex polarFront(const imageComplex &img, bool logFlag = false);

imageReal<int> circularMask(const long long shape[2], double radius);

imageComplex bind(const imageReal<double> &abs, const imageReal<double> &angle);

#endif // ALIGNMENT_IMAGE2D_H
