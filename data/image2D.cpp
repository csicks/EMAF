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

#include "image2D.h"

imageComplex::imageComplex() {
    shape[0] = 0;
    shape[1] = 0;
    data = arrayComplex();
}

imageComplex::imageComplex(const int s[2]) {
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayComplex(static_cast<size_t>(shape[0]) * shape[1]);
}

imageComplex::imageComplex(const imageComplex &img) {
    shape[0] = img.shape[0];
    shape[1] = img.shape[1];
    data = arrayComplex(img.data);
}

imageComplex::imageComplex(const imageComplex *img) : imageComplex(*img) {
}

imageComplex::imageComplex(const arrayComplex &d, const int s[2]) {
    if (d.length != static_cast<size_t>(s[0]) * s[1])
        throw baseException("Error: Data and shape are not matched!");
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayComplex(d);
}

imageComplex::imageComplex(const double d[][2], const int s[2]) {
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayComplex(d, static_cast<size_t>(s[0]) * s[1]);
}

double *imageComplex::get(int i, int j) const {
    return data[static_cast<size_t>(i) * shape[1] + j];
}

imageComplex imageComplex::get(int i, int j, int length, int width) const {
    if (i + length > shape[0] or j + width > shape[1])
        throw baseException("Error: Index out of range!");
    int s[2] = {length, width};
    imageComplex r = imageComplex(s);
    for (int m = 0; m < length; ++m)
        for (int n = 0; n < width; ++n) {
            double *value = get(m + i, n + j);
            r.set(m, n, value);
            delete[] value;
        }
    return r;
}

void imageComplex::set(int i, int j, const double value[2]) {
    data.set(static_cast<size_t>(i) * shape[1] + j, value);
}

void imageComplex::set(int i, int j, const imageComplex &img) {
    if (i + img.shape[0] > shape[0] or j + img.shape[1] > shape[1])
        throw baseException("Error: Index out of range!");
    for (int m = 0; m < img.shape[0]; ++m)
        for (int n = 0; n < img.shape[1]; ++n) {
            double *value = img.get(m, n);
            set(m + i, n + j, value);
            delete[] value;
        }
}

double *imageComplex::mean() const {
    return data.mean();
}

imageComplex imageComplex::copy() const {
    imageComplex img = imageComplex(this);
    return img;
}

imageReal<double> imageComplex::abs() const {
    imageReal<double> img = imageReal<double>(shape);
    img.data = data.abs();
    return img;
}

imageReal<double> imageComplex::real() const {
    imageReal<double> img = imageReal<double>(shape);
    img.data = data.real();
    return img;
}

imageReal<double> imageComplex::imag() const {
    imageReal<double> img = imageReal<double>(shape);
    img.data = data.imag();
    return img;
}

imageReal<double> imageComplex::angle() const {
    imageReal<double> img = imageReal<double>(shape);
    img.data = data.angle();
    return img;
}

arrayComplex imageComplex::operator[](int index) const {
    arrayComplex ary = arrayComplex(shape[1]);
    memcpy(ary.data, data.data + 2 * static_cast<size_t>(index) * shape[1], 2 * sizeof(double) * shape[1]);
    return ary;
}

imageComplex &imageComplex::operator=(const imageComplex &img) {
    if (this == &img) {
        return *this;
    }
    shape[0] = img.shape[0];
    shape[1] = img.shape[1];
    data = img.data;
    return *this;
}

imageComplex imageComplex::operator+(double value) const {
    imageComplex img = imageComplex(shape);
    img.data = data + value;
    return img;
}

imageComplex imageComplex::operator-(double value) const {
    imageComplex img = imageComplex(shape);
    img.data = data - value;
    return img;
}

imageComplex imageComplex::operator*(double value) const {
    imageComplex img = imageComplex(shape);
    img.data = data * value;
    return img;
}

imageComplex imageComplex::operator/(double value) const {
    imageComplex img = imageComplex(shape);
    img.data = data / value;
    return img;
}

imageComplex imageComplex::operator+(const double value[2]) const {
    imageComplex img = imageComplex(shape);
    img.data = data + value;
    return img;
}

imageComplex imageComplex::operator-(const double value[2]) const {
    imageComplex img = imageComplex(shape);
    img.data = data - value;
    return img;
}

imageComplex imageComplex::operator*(const double value[2]) const {
    imageComplex img = imageComplex(shape);
    img.data = data * value;
    return img;
}

imageComplex imageComplex::operator/(const double value[2]) const {
    imageComplex img = imageComplex(shape);
    img.data = data / value;
    return img;
}

void imageComplex::print() const {
    std::cout << 'c' << '\n';
    for (int i = 0; i < shape[0]; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            double *r = get(i, j);
            std::cout << r[0];
            if (r[1] >= 0)
                std::cout << '+' << r[1] << 'i' << ' ';
            else
                std::cout << r[1] << 'i' << ' ';
            delete[] r;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

imageComplex imageComplex::conj() const {
    imageComplex img = imageComplex(shape);
    img.data = data.conj();
    return img;
}

imageComplex imageComplex::operator+(const imageComplex &img) const {
    if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
        throw baseException("Error: Two images are not of the same size!");
    imageComplex imgOut = imageComplex(shape);
    imgOut.data = data + img.data;
    return imgOut;
}

imageComplex imageComplex::operator-(const imageComplex &img) const {
    if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
        throw baseException("Error: Two images are not of the same size!");
    imageComplex imgOut = imageComplex(shape);
    imgOut.data = data - img.data;
    return imgOut;
}

imageComplex imageComplex::operator*(const imageComplex &img) const {
    if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
        throw baseException("Error: Two images are not of the same size!");
    imageComplex imgOut = imageComplex(shape);
    imgOut.data = data * img.data;
    return imgOut;
}

imageComplex imageComplex::operator/(const imageComplex &img) const {
    if (shape[0] != img.shape[0] or shape[1] != img.shape[1])
        throw baseException("Error: Two images are not of the same size!");
    imageComplex imgOut = imageComplex(shape);
    imgOut.data = data / img.data;
    return imgOut;
}

imageComplex imageComplex::flipX() const {
    imageComplex img = imageComplex(shape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j) {
            double *value = get(static_cast<int>(shape[0]) - 1 - i, j);
            img.set(i, j, value);
            delete[] value;
        }
    return img;
}

imageComplex imageComplex::flipY() const {
    imageComplex img = imageComplex(shape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j) {
            double *value = get(i, static_cast<int>(shape[1]) - 1 - j);
            img.set(i, j, value);
            delete[] value;
        }
    return img;
}

imageComplex imageComplex::flipAll() const {
    imageComplex img = imageComplex(shape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j) {
            double *value = get(static_cast<int>(shape[0]) - 1 - i, static_cast<int>(shape[1]) - 1 - j);
            img.set(i, j, value);
            delete[] value;
        }
    return img;
}

arrayComplex imageComplex::row(int index) const {
    arrayComplex ary = arrayComplex(shape[1]);
    for (int i = 0; i < shape[1]; ++i) {
        double *value = get(index, i);
        ary.set(i, value);
        delete[] value;
    }
    return ary;
}

arrayComplex imageComplex::column(int index) const {
    arrayComplex ary = arrayComplex(shape[0]);
    for (int i = 0; i < shape[1]; ++i) {
        double *value = get(i, index);
        ary.set(i, value);
        delete[] value;
    }
    return ary;
}

double *interBiLinear(const imageComplex &img, double x, double y) {
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    arrayComplex neighbour = arrayComplex();
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
    auto *result = new double[2];

    if (inf.sum() == 0) {
        result[0] = 0;
        result[1] = 0;
    } else if (inf.sum() == 1) {
        for (int i = 0; i < inf.length; ++i)
            if (inf[i] == 1) {
                result = neighbour[i];
                break;
            }
    } else {
        if (inf.equalList({1, 1, 0, 0})) {
            result[0] = neighbour[0][0] * (y - yl) + neighbour[1][0] * (yr - y);
            result[1] = neighbour[0][1] * (y - yl) + neighbour[1][1] * (yr - y);
        } else if (inf.equalList({1, 0, 1, 0})) {
            result[0] = neighbour[0][0] * (x - xl) + neighbour[2][0] * (xr - x);
            result[1] = neighbour[0][1] * (x - xl) + neighbour[2][1] * (xr - x);
        } else if (inf.equalList({0, 1, 0, 1})) {
            result[0] = neighbour[1][0] * (x - xl) + neighbour[3][0] * (xr - x);
            result[1] = neighbour[1][1] * (x - xl) + neighbour[3][1] * (xr - x);
        } else if (inf.equalList({0, 0, 1, 1})) {
            result[0] = neighbour[2][0] * (y - yl) + neighbour[3][0] * (yr - y);
            result[1] = neighbour[2][1] * (y - yl) + neighbour[3][1] * (yr - y);
        } else if (inf.equalList({1, 1, 1, 0})) {
            double temp0 = neighbour[0][0] * (y - yl) + neighbour[1][0] * (yr - y);
            double temp1 = neighbour[0][1] * (y - yl) + neighbour[1][1] * (yr - y);
            result[0] = temp0 * (x - xl) + neighbour[2][0] * (xr - x);
            result[1] = temp1 * (x - xl) + neighbour[2][1] * (xr - x);
        } else if (inf.equalList({1, 1, 0, 1})) {
            double temp0 = neighbour[0][0] * (y - yl) + neighbour[1][0] * (yr - y);
            double temp1 = neighbour[0][1] * (y - yl) + neighbour[1][1] * (yr - y);
            result[0] = temp0 * (x - xl) + neighbour[3][0] * (xr - x);
            result[1] = temp1 * (x - xl) + neighbour[3][1] * (xr - x);
        } else if (inf.equalList({1, 0, 1, 1})) {
            double temp0 = neighbour[2][0] * (y - yl) + neighbour[3][0] * (yr - y);
            double temp1 = neighbour[2][1] * (y - yl) + neighbour[3][1] * (yr - y);
            result[0] = neighbour[0][0] * (x - xl) + temp0 * (xr - x);
            result[1] = neighbour[0][1] * (x - xl) + temp1 * (xr - x);
        } else if (inf.equalList({0, 1, 1, 1})) {
            double temp0 = neighbour[2][0] * (y - yl) + neighbour[3][0] * (yr - y);
            double temp1 = neighbour[2][1] * (y - yl) + neighbour[3][1] * (yr - y);
            result[0] = neighbour[1][0] * (x - xl) + temp0 * (xr - x);
            result[1] = neighbour[1][1] * (x - xl) + temp1 * (xr - x);
        } else {
            double temp10 = neighbour[0][0] * (y - yl) + neighbour[1][0] * (yr - y);
            double temp11 = neighbour[0][1] * (y - yl) + neighbour[1][1] * (yr - y);
            double temp20 = neighbour[2][0] * (y - yl) + neighbour[3][0] * (yr - y);
            double temp21 = neighbour[2][1] * (y - yl) + neighbour[3][1] * (yr - y);
            result[0] = temp10 * (x - xl) + temp20 * (xr - x);
            result[1] = temp11 * (x - xl) + temp21 * (xr - x);
        }
    }
    return result;
}

double *interAvg(const imageComplex &img, double x, double y) {
    int xl = std::floor(x);
    int xr = std::ceil(x);
    int yl = std::floor(y);
    int yr = std::ceil(y);
    arrayComplex neighbour = arrayComplex();
    if (0 <= xl and xl < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xl, yl));
    if (0 <= xl and xl < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xl, yr));
    if (0 <= xr and xr < img.shape[0] and 0 <= yl and yl < img.shape[1])
        neighbour.append(img.get(xr, yl));
    if (0 <= xr and xr < img.shape[0] and 0 <= yr and yr < img.shape[1])
        neighbour.append(img.get(xr, yr));
    auto *result = new double[2];
    if (neighbour.empty()) {
        result[0] = 0;
        result[1] = 0;
    } else {
        double *s = neighbour.sum();
        result[0] = s[0] / static_cast<double>(neighbour.length);
        result[1] = s[1] / static_cast<double>(neighbour.length);
        delete[] s;
    }
    return result;
}

imageComplex warp(const imageComplex &img, double matrix[3][3]) {
    imageComplex imgNew = imageComplex(img.shape);
    double abs = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double matrix_inv[3][3] = {{matrix[1][1] / abs,  -matrix[0][1] / abs,
                                                                         (matrix[0][1] * matrix[1][2] -
                                                                          matrix[1][1] * matrix[0][2]) / abs},
                               {-matrix[1][0] / abs, matrix[0][0] / abs, (matrix[1][0] * matrix[0][2] -
                                                                          matrix[0][0] * matrix[1][2]) / abs},
                               {0,                   0,                  1}};
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            double x_in = matrix_inv[0][0] * i + matrix_inv[0][1] * j + matrix_inv[0][2];
            double y_in = matrix_inv[1][0] * i + matrix_inv[1][1] * j + matrix_inv[1][2];
            imgNew.set(i, j, interAvg(img, x_in, y_in));
        }
    return imgNew;
}

imageComplex shift(const imageComplex &img, double x, double y) {
    double matrix[3][3] = {{1, 0, x},
                           {0, 1, y},
                           {0, 0, 1}};
    imageComplex imgNew = warp(img, matrix);
    return imgNew;
}

imageComplex rotate(const imageComplex &img, double angle) {
    double angleR = -angle / 360 * 2 * M_PI;
    double rows = img.shape[0] / 2.0;
    double cols = img.shape[1] / 2.0;
    double matrix[3][3] = {{std::cos(angleR),  std::sin(angleR),
                                                                 rows - rows * std::cos(angleR) -
                                                                 cols * std::sin(angleR)},
                           {-std::sin(angleR), std::cos(angleR), cols - cols * std::cos(angleR) +
                                                                 rows * std::sin(angleR)},
                           {0,                 0,                1}};
    imageComplex imgNew = warp(img, matrix);
    return imgNew;
}

imageComplex hanning(const imageComplex &img) {
    arrayReal<double> row = arrayReal<double>(img.shape[1]);
    arrayReal<double> col = arrayReal<double>(img.shape[0]);
    for (int i = 0; i < img.shape[1]; ++i)
        row[i] = 0.5 - 0.5 * std::cos(2 * M_PI * i / static_cast<int>(img.shape[1] - 1));
    for (int i = 0; i < img.shape[0]; ++i)
        col[i] = 0.5 - 0.5 * std::cos(2 * M_PI * i / static_cast<int>(img.shape[1] - 1));
    imageReal<double> imgRow = extend2Image(row, img.shape[0], 0);
    imageReal<double> imgCol = extend2Image(col, img.shape[1], 1);
    imageReal<double> mask = imgRow * imgCol;
    return img * mask;
}

imageComplex hamming(const imageComplex &img) {
    arrayReal<double> row = arrayReal<double>(img.shape[1]);
    arrayReal<double> col = arrayReal<double>(img.shape[0]);
    for (int i = 0; i < img.shape[1]; ++i)
        row[i] = 0.54 - 0.46 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    for (int i = 0; i < img.shape[0]; ++i)
        col[i] = 0.54 - 0.46 * std::cos(2 * M_PI * i / (img.shape[1] - 1));
    imageReal<double> imgRow = extend2Image(row, img.shape[0], 0);
    imageReal<double> imgCol = extend2Image(col, img.shape[1], 1);
    imageReal<double> mask = imgRow * imgCol;
    return img * mask;
}

imageComplex butterworthLow(const imageComplex &img, const int center[2], double d0, int n) {
    imageComplex r = imageComplex(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            double dn = std::sqrt((i - center[0]) * (i - center[0]) + (j - center[1]) * (j - center[1]));
            double rate = 1 + std::pow(dn / d0, n);
            double *value = img.get(i, j);
            double valueN[2] = {value[0] / rate, value[1] / rate};
            delete[] value;
            r.set(i, j, valueN);
        }
    return r;
}

imageComplex butterworthHigh(const imageComplex &img, const int center[2], double d0, int n) {
    imageComplex r = imageComplex(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            double dn = std::sqrt((i - center[0]) * (i - center[0]) + (j - center[1]) * (j - center[1]));
            if (dn == 0) {
                double valueN[2] = {0, 0};
                r.set(i, j, valueN);
            } else {
                double rate = 1 + std::pow(d0 / dn, n);
                double *value = img.get(i, j);
                double valueN[2] = {value[0] / rate, value[1] / rate};
                delete[] value;
                r.set(i, j, valueN);
            }
        }
    return r;
}

imageComplex polarBack(const imageComplex &img, bool logFlag) {
    double centerX = img.shape[0] / 2.0;
    double centerY = img.shape[1] / 2.0;
    imageComplex imgPolar = imageComplex(img.shape);
    double maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            double angle = static_cast<double>(i) / img.shape[0] * 360;
            double radius;
            if (logFlag)
                radius = std::exp(static_cast<double>(j) / img.shape[1] * std::log(maxRadius));
            else
                radius = static_cast<double>(j) / img.shape[1] * maxRadius;
            //            double x = std::round(-radius * std::sin(angle / 360 * M_PI * 2) + centerX);
            //            double y = std::round(-radius * std::cos(angle / 360 * M_PI * 2) + centerY);
            double x = -radius * std::sin(angle / 360 * M_PI * 2) + centerX;
            double y = -radius * std::cos(angle / 360 * M_PI * 2) + centerY;
            double *value = interAvg(img, x, y);
            imgPolar.set(i, j, value);
            delete[] value;
        }
    return imgPolar;
}

arrayComplex findNeighbour(const imageComplex &img, int x, int y) {
    int xl = x - 1;
    int xr = x + 1;
    int yl = y - 1;
    int yr = y + 1;
    arrayComplex neighbour = arrayComplex();
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

imageComplex valueInter(const imageComplex &img, const double value[2]) {
    imageComplex imgNew = img.copy();
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
            if (std::sqrt(std::pow(img.get(i, j)[0] - value[0], 2) + std::pow(img.get(i, j)[1] - value[1], 2)) <
                MIN_D) {
                double *v = findNeighbour(img, i, j).mean();
                imgNew.set(i, j, v);
                delete[] v;
            }
    return imgNew;
}

imageComplex polarFront(const imageComplex &img, bool logFlag) {
    double centerX = img.shape[0] / 2.0;
    double centerY = img.shape[1] / 2.0;
    imageComplex imgPolar = imageComplex(img.shape);
    double maxRadius = std::sqrt(std::pow(centerX, 2) + std::pow(centerY, 2));
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            double radius = std::sqrt(std::pow(i - centerX, 2) + std::pow(j - centerY, 2));
            double yD = centerX - i;
            double xD = j - centerY;
            double angle;
            if (xD == 0 and yD > 0)
                angle = 90;
            else if (xD == 0 and yD < 0)
                angle = 270;
            else if (xD > 0 and yD >= 0)
                angle = std::atan(yD / xD) / 2 / M_PI * 360;
            else if (xD < 0)
                angle = std::atan(yD / xD) / 2 / M_PI * 360 + 180;
            else if (xD > 0 and yD < 0)
                angle = std::atan(yD / xD) / 2 / M_PI * 360 + 360;
            else
                angle = 0;
            double radiusIndex;
            if (logFlag)
                radiusIndex = radius == 0 ? 0 : std::round(std::log(radius) / std::log(maxRadius) * (img.shape[1] - 1));
            else
                radiusIndex = std::round(radius / maxRadius * (img.shape[1] - 1));
            double angIndex = -std::round(angle / 360 * (img.shape[0] - 1));
            double *value = img.get(i, j);
            imgPolar.set(static_cast<int>(angIndex), static_cast<int>(radiusIndex), value);
            delete[] value;
        }
    double value[2] = {0, 0};
    for (int i = 0; i < 3; ++i)
        imgPolar = valueInter(imgPolar, value);
    return imgPolar;
}

imageReal<int> circularMask(const int s[2], double radius) {
    double centerX = s[0] / 2.0;
    double centerY = s[1] / 2.0;
    imageReal<int> r = imageReal<int>(s);
    for (int i = 0; i < s[0]; ++i)
        for (int j = 0; j < s[1]; ++j) {
            if ((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY) <= radius * radius)
                r.set(i, j, 1);
            else
                r.set(i, j, 0);
        }
    return r;
}

imageComplex bind(const imageReal<double> &abs, const imageReal<double> &angle) {
    imageComplex r = imageComplex(abs.shape);
    r.data = bind(abs.data, angle.data);
    return r;
}