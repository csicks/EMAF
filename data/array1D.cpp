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

#include "array1D.h"

arrayComplex::arrayComplex()
{
    data = new double[2 * maxLength];
    length = 0;
    maxLength = maxLengthInit;
}

arrayComplex::arrayComplex(const arrayComplex &ary)
{
    data = new double[2 * ary.maxLength];
    length = ary.length;
    maxLength = ary.maxLength;
    memcpy(data, ary.data, 2 * length * sizeof(double));
}

arrayComplex::arrayComplex(const arrayComplex *ary) : arrayComplex(*ary)
{
}

arrayComplex::arrayComplex(const double d[][2], int l)
{
    maxLength = l > maxLengthInit ? l : maxLengthInit;
    data = new double[2 * maxLength];
    length = l;
    for (int i = 0; i < l; ++i)
    {
        data[2 * i] = d[i][0];
        data[2 * i + 1] = d[i][1];
    }
}

arrayComplex::arrayComplex(int l)
{
    maxLength = l > maxLengthInit ? l : maxLengthInit;
    data = new double[2 * maxLength];
    for (int i = 0; i < 2 * l; ++i)
        data[i] = 0;
    length = l;
}

arrayComplex::~arrayComplex()
{
    if (data)
    {
        delete[] data;
        data = nullptr;
    }
}

void arrayComplex::append(const double value[2])
{
    if (length + 1 > maxLength)
    {
        auto *newData = new double[2 * length + 2];
        memcpy(newData, data, 2 * length * sizeof(double));
        newData[2 * length] = value[0];
        newData[2 * length + 1] = value[1];
        delete[] data;
        data = newData;
        length++;
        maxLength++;
    }
    else
    {
        data[2 * length] = value[0];
        data[2 * length + 1] = value[1];
        length++;
    }
}

void arrayComplex::append(double value)
{
    if (length + 1 > maxLength)
    {
        auto *newData = new double[2 * length + 2];
        memcpy(newData, data, 2 * length * sizeof(double));
        newData[2 * length] = value;
        newData[2 * length + 1] = 0;
        delete[] data;
        data = newData;
        length++;
        maxLength++;
    }
    else
    {
        data[2 * length] = value;
        data[2 * length + 1] = 0;
        length++;
    }
}

void arrayComplex::remove()
{
    if (length > maxLengthInit)
    {
        maxLength = length - 1 > maxLengthInit ? length - 1 : maxLengthInit;
        auto *newData = new double[2 * maxLength];
        memcpy(newData, data, 2 * (length - 1) * sizeof(double));
        delete[] data;
        data = newData;
        length--;
    }
    else if (length == 0)
        throw baseException("Error: Array is empty!");
    else
        length--;
}

double *arrayComplex::get(int index) const
{
    auto *r = new double[2];
    r[0] = data[2 * index];
    r[1] = data[2 * index + 1];
    return r;
}

arrayComplex arrayComplex::get(int index, int l) const
{
    if (index + l > length)
        throw baseException("Error: Index out of range!");
    arrayComplex ary = arrayComplex(l);
    memcpy(ary.data, data + 2 * index, 2 * l * sizeof(double));
    return ary;
}

void arrayComplex::set(int index, const arrayComplex &ary)
{
    if (index + ary.length > length)
        throw baseException("Error: Index out of range!");
    memcpy(data + 2 * index, ary.data, 2 * ary.length * sizeof(double));
}

void arrayComplex::set(int index, const double value[2])
{
    if (index > length - 1)
        throw baseException("Error: Index out of range!");
    data[2 * index] = value[0];
    data[2 * index + 1] = value[1];
}

arrayComplex arrayComplex::copy() const
{
    arrayComplex ary = arrayComplex(this);
    return ary;
}

bool arrayComplex::empty() const
{
    return length == 0;
}

arrayReal<double> arrayComplex::abs() const
{
    arrayReal<double> ary = arrayReal<double>(length);
    for (int i = 0; i < length; ++i)
        ary[i] = std::sqrt(data[2 * i] * data[2 * i] + data[2 * i + 1] * data[2 * i + 1]);
    return ary;
}

arrayReal<double> arrayComplex::angle() const
{
    arrayReal<double> ary = arrayReal<double>(length);
    for (int i = 0; i < length; ++i)
        ary[i] = std::atan2(data[2 * i + 1], data[2 * i]);
    return ary;
}

arrayReal<double> arrayComplex::real() const
{
    arrayReal<double> ary = arrayReal<double>(length);
    for (int i = 0; i < length; ++i)
        ary[i] = data[2 * i];
    return ary;
}

arrayReal<double> arrayComplex::imag() const
{
    arrayReal<double> ary = arrayReal<double>(length);
    for (int i = 0; i < length; ++i)
        ary[i] = data[2 * i + 1];
    return ary;
}

double *arrayComplex::sum() const
{
    auto *s = new double[2];
    for (int i = 0; i < length; ++i)
    {
        s[0] += data[2 * i];
        s[1] += data[2 * i + 1];
    }
    return s;
}

double *arrayComplex::mean() const
{
    auto *r = new double[2];
    for (int i = 0; i < length; ++i)
    {
        r[0] += data[2 * i] / length;
        r[1] += data[2 * i + 1] / length;
    }
    return r;
}

double *arrayComplex::operator[](int index) const
{
    auto *r = new double[2];
    r[0] = data[2 * index];
    r[1] = data[2 * index + 1];
    return r;
}

arrayComplex &arrayComplex::operator=(const arrayComplex &ary)
{
    if (this == &ary)
    {
        return *this;
    }
    if (length == ary.length)
    {
        memcpy(data, ary.data, 2 * length * sizeof(double));
        return *this;
    }
    else
    {
        delete[] data;
        data = new double[2 * ary.maxLength];
        length = ary.length;
        maxLength = ary.maxLength;
        memcpy(data, ary.data, 2 * length * sizeof(double));
        return *this;
    }
}

arrayComplex arrayComplex::conj() const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < length; ++i)
        ary.data[2 * i + 1] *= -1;
    return ary;
}

void arrayComplex::print() const
{
    std::cout << 'c' << '\t';
    for (int i = 0; i < length; ++i)
    {
        std::cout << data[2 * i];
        if (data[2 * i + 1] >= 0)
            std::cout << '+' << data[2 * i + 1] << 'i' << ' ';
        else
            std::cout << data[2 * i + 1] << 'i' << ' ';
    }
    std::cout << std::endl;
}

arrayComplex arrayComplex::operator+(double value) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < length; ++i)
        ary.data[2 * i] += value;
    return ary;
}

arrayComplex arrayComplex::operator-(double value) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < length; ++i)
        ary.data[2 * i] -= value;
    return ary;
}

arrayComplex arrayComplex::operator*(double value) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < 2 * length; ++i)
        ary.data[i] *= value;
    return ary;
}

arrayComplex arrayComplex::operator/(double value) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < 2 * length; ++i)
        ary.data[i] = value == 0 ? INFINITY : ary.data[i] / value;
    return ary;
}

arrayComplex arrayComplex::operator+(const double value[2]) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < length; ++i)
    {
        ary.data[2 * i] += value[0];
        ary.data[2 * i + 1] += value[1];
    }
    return ary;
}

arrayComplex arrayComplex::operator-(const double value[2]) const
{
    arrayComplex ary = arrayComplex(this);
    for (int i = 0; i < length; ++i)
    {
        ary.data[2 * i] -= value[0];
        ary.data[2 * i + 1] -= value[1];
    }
    return ary;
}

arrayComplex arrayComplex::operator*(const double value[2]) const
{
    arrayComplex ary = arrayComplex(length);
    for (int i = 0; i < length; ++i)
    {
        ary.data[2 * i] = data[2 * i] * value[0] - data[2 * i + 1] * value[1];
        ary.data[2 * i + 1] = data[2 * i] * value[1] + data[2 * i + 1] * value[0];
    }
    return ary;
}

arrayComplex arrayComplex::operator/(const double value[2]) const
{
    arrayComplex ary = arrayComplex(this);
    double c[2] = {value[0], -value[1]};
    ary = ary * c;
    ary = ary / (value[0] * value[0] + value[1] * value[1]);
    return ary;
}

arrayReal<int> arrayComplex::isNotInf() const
{
    arrayReal<int> index = arrayReal<int>(length);
    for (int i = 0; i < length; ++i)
        index[i] = data[2 * i] == INFINITY ? 0 : 1;
    return index;
}

arrayComplex arrayComplex::operator+(const arrayComplex &ary) const
{
    if (ary.length != length)
        throw baseException("Error: Two arrays are not of the same shape!");
    arrayComplex aryOut = arrayComplex(this);
    for (int i = 0; i < 2 * length; ++i)
        aryOut.data[i] += ary.data[i];
    return aryOut;
}

arrayComplex arrayComplex::operator-(const arrayComplex &ary) const
{
    if (ary.length != length)
        throw baseException("Error: Two arrays are not of the same shape!");
    arrayComplex aryOut = arrayComplex(this);
    for (int i = 0; i < 2 * length; ++i)
        aryOut.data[i] -= ary.data[i];
    return aryOut;
}

arrayComplex arrayComplex::operator*(const arrayComplex &ary) const
{
    if (ary.length != length)
        throw baseException("Error: Two arrays are not of the same shape!");
    arrayComplex aryOut = arrayComplex(length);
    for (int i = 0; i < length; ++i)
    {
        aryOut.data[2 * i] = data[2 * i] * ary.data[2 * i] - data[2 * i + 1] * ary.data[2 * i + 1];
        aryOut.data[2 * i + 1] = data[2 * i] * ary.data[2 * i + 1] + data[2 * i + 1] * ary.data[2 * i];
    }
    return aryOut;
}

arrayComplex arrayComplex::operator/(const arrayComplex &ary) const
{
    if (ary.length != length)
        throw baseException("Error: Two arrays are not of the same shape!");
    arrayComplex aryOut = arrayComplex(length);
    for (int i = 0; i < length; ++i)
    {
        double value = (ary.data[2 * i] * ary.data[2 * i] + ary.data[2 * i + 1] * ary.data[2 * i + 1]);
        if (value == 0)
        {
            aryOut.data[2 * i] = INFINITY;
            aryOut.data[2 * i + 1] = 0;
        }
        else
        {
            aryOut.data[2 * i] = (data[2 * i] * ary.data[2 * i] + data[2 * i + 1] * ary.data[2 * i + 1]) / value;
            aryOut.data[2 * i + 1] = (-data[2 * i] * ary.data[2 * i + 1] + data[2 * i + 1] * ary.data[2 * i]) / value;
        }
    }
    return aryOut;
}

arrayComplex bind(const arrayReal<double> &abs, const arrayReal<double> &angle)
{
    arrayComplex r = arrayComplex(abs.length);
    for (int i = 0; i < r.length; ++i)
    {
        r.data[2 * i] = abs[i] * std::cos(angle[i]);
        r.data[2 * i + 1] = abs[i] * std::sin(angle[i]);
    }
    return r;
}