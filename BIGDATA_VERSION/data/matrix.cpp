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

#include "matrix.h"

matrix::matrix()
{
    data = arrayReal<double>();
}

matrix::matrix(const matrix &matrix)
{
    data = matrix.data;
}

matrix::matrix(const int s[2])
{
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayReal<double>(s[0] * s[1]);
}

matrix::matrix(const std::initializer_list<double> &initList, const int s[2])
{
    if (initList.size() != s[0] * s[1])
        throw baseException("Error: Data and shape are not matched!");
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayReal<double>(initList);
}

matrix::matrix(const arrayReal<double> &d, const int *s)
{
    if (d.length != s[0] * s[1])
        throw baseException("Error: Data and shape are not matched!");
    shape[0] = s[0];
    shape[1] = s[1];
    data = arrayReal<double>(d);
}

matrix::matrix(const double *d, const int s[2])
{
    data = arrayReal<double>(d, s[0] * s[1]);
    shape[0] = s[0];
    shape[1] = s[1];
}

matrix &matrix::operator=(const matrix &matrix)
{
    if (this == &matrix)
    {
        return *this;
    }
    data = matrix.data;
    shape[0] = matrix.shape[0];
    shape[1] = matrix.shape[1];
    return *this;
}

double matrix::get(int i, int j) const
{
    return data[i * shape[1] + j];
}

void matrix::set(int i, int j, double value)
{
    data[i * shape[1] + j] = value;
}

arrayReal<double> matrix::operator[](int index) const
{
    arrayReal<double> ary = arrayReal<double>(shape[1]);
    memcpy(ary.data, data.data + index * shape[1], sizeof(double) * shape[1]);
    return ary;
}

arrayReal<double> matrix::row(int index) const
{
    arrayReal<double> ary = arrayReal<double>(shape[1]);
    memcpy(ary.data, data.data + index * shape[1], sizeof(double) * shape[1]);
    return ary;
}

arrayReal<double> matrix::column(int index) const
{
    arrayReal<double> ary = arrayReal<double>(shape[0]);
    for (int i = 0; i < shape[0]; ++i)
        ary[i] = get(i, index);
    return ary;
}

matrix matrix::removeRow(int index) const
{
    int newShape[2] = {shape[0] - 1, shape[1]};
    matrix temp = matrix(newShape);
    memcpy(temp.data.data, data.data, index * shape[1] * sizeof(double));
    memcpy(temp.data.data + index * shape[1], data.data + (index + 1) * shape[1], (shape[0] - index - 1) * shape[1] * sizeof(double));
    return temp;
}

matrix matrix::removeColumn(int index) const
{
    int newShape[2] = {shape[0], shape[1] - 1};
    matrix temp = matrix(newShape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j)
        {
            if (j < index)
            {
                temp.set(i, j, get(i, j));
            }
            else if (j > index)
            {
                temp.set(i, j - 1, get(i, j));
            }
        }
    return temp;
}

matrix matrix::operator+(const matrix &matrixIn) const
{
    if (shape[0] != matrixIn.shape[0] || shape[1] != matrixIn.shape[1])
        throw baseException("Error: Two matrices are not of the same size!");
    matrix matrixOut = matrix(shape);
    matrixOut.data = data + matrixIn.data;
    return matrixOut;
}

matrix matrix::operator-(const matrix &matrixIn) const
{
    if (shape[0] != matrixIn.shape[0] || shape[1] != matrixIn.shape[1])
        throw baseException("Error: Two matrices are not of the same size!");
    matrix matrixOut = matrix(shape);
    matrixOut.data = data - matrixIn.data;
    return matrixOut;
}

matrix matrix::operator*(const matrix &matrixIn) const
{
    matrix matrixOut = matrix(shape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j)
        {
            arrayReal<double> r = row(i);
            arrayReal<double> c = matrixIn.column(j);
            double value = (r * c).sum();
            matrixOut.set(i, j, value);
        }
    return matrixOut;
}

double matrix::determinant() const
{
    if (shape[0] != shape[1])
        throw baseException("Error: Only square matrices have determinants!");
    if (shape[0] == 2)
        return get(0, 0) * get(1, 1) - get(0, 1) * get(1, 0);
    else
    {
        double s = 0;
        matrix temp = removeRow(0);
        for (int i = 0; i < shape[1]; ++i)
        {
            double flag = std::pow(-1, i) * get(0, i);
            matrix remain = temp.removeColumn(i);
            s += flag * remain.determinant();
        }
        return s;
    }
}

matrix matrix::transpose() const
{
    int s[2] = {shape[1], shape[0]};
    matrix matrixOut = matrix(s);
    for (int i = 0; i < s[0]; ++i)
        for (int j = 0; j < s[1]; ++j)
        {
            matrixOut.set(i, j, get(j, i));
        }
    return matrixOut;
}

matrix matrix::adjoint() const
{
    matrix matrixOut = matrix(shape);
    for (int i = 0; i < shape[0]; ++i)
        for (int j = 0; j < shape[1]; ++j)
        {
            matrix temp = removeRow(i);
            temp = temp.removeColumn(j);
            double value = temp.determinant() * std::pow(-1, i + j);
            matrixOut.set(i, j, value);
        }
    matrixOut = matrixOut.transpose();
    return matrixOut;
}

matrix matrix::inverse() const
{
    double value = determinant();
    if (value == 0)
        throw baseException("Error: Inverse matrix does not exist!");
    matrix matrixOut = adjoint();
    matrixOut = matrixOut / value;
    return matrixOut;
}

void matrix::print() const
{
    for (int i = 0; i < shape[0]; ++i)
    {
        for (int j = 0; j < shape[1]; ++j)
            std::cout << get(i, j) << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

matrixTransform::matrixTransform(double angle, int x, int y)
{
    shape[0] = shape[1] = 3;
    data = arrayReal<double>(shape[0] * shape[1]);
    angle = angle / 360 * 2 * M_PI;
    data[0] = std::cos(angle);
    data[1] = std::sin(angle);
    data[2] = x;
    data[3] = -std::sin(angle);
    data[4] = std::cos(angle);
    data[5] = y;
    data[6] = data[7] = 0;
    data[8] = 1;
}

matrixTransform::matrixTransform(const std::initializer_list<double> &initList)
{
    if (initList.size() != 9)
        throw baseException("Error: Data and shape are not matched!");
    shape[0] = shape[1] = 3;
    data = arrayReal<double>(initList);
}

matrixTransform::matrixTransform(const matrix &m)
{
    if (m.shape[0] != 3 || m.shape[1] != 3)
        throw baseException("Error: Shape of transform matrix should be 3*3!");
    shape[0] = m.shape[0];
    shape[1] = m.shape[1];
    data = m.data;
}

matrixTransform &matrixTransform::operator=(const matrixTransform &m)
{
    if (this == &m)
    {
        return *this;
    }
    data = m.data;
    return *this;
}