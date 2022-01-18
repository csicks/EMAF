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

#ifndef ALIGNMENT_MATRIX_H
#define ALIGNMENT_MATRIX_H

#include "array1D.h"

class matrix
{
public:
    int shape[2]{};
    arrayReal<double> data;

public:
    matrix();

    matrix(const matrix &matrix);

    explicit matrix(const int s[2]);

    matrix(const std::initializer_list<double> &initList, const int s[2]);

    matrix(const arrayReal<double> &d, const int *s);

    matrix(const double *d, const int s[2]);

    matrix &operator=(const matrix &matrix);

    double get(int i, int j) const;

    void set(int i, int j, double value);

    arrayReal<double> operator[](int index) const;

    arrayReal<double> row(int index) const;

    arrayReal<double> column(int index) const;

    matrix removeRow(int index) const;

    matrix removeColumn(int index) const;

    matrix operator+(const matrix &matrixIn) const;

    matrix operator-(const matrix &matrixIn) const;

    matrix operator*(const matrix &matrixIn) const;

    template <typename T>
    matrix operator+(T value) const
    {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data + value;
        return matrixOut;
    }

    template <typename T>
    matrix operator-(T value) const
    {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data - value;
        return matrixOut;
    }

    template <typename T>
    matrix operator*(T value) const
    {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data * value;
        return matrixOut;
    }

    template <typename T>
    matrix operator/(T value) const
    {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data / value;
        return matrixOut;
    }

    double determinant() const;

    matrix transpose() const;

    matrix adjoint() const;

    matrix inverse() const;

    void print() const;
};

/// matrix used for transformation estimation, shape should be 3*3
class matrixTransform : public matrix
{

public:
    using matrix::matrix;

    matrixTransform(double angle, int x, int y);

    matrixTransform(const std::initializer_list<double> &initList);

    explicit matrixTransform(const matrix &m);

    matrixTransform &operator=(const matrixTransform &m);
};

#endif // ALIGNMENT_MATRIX_H