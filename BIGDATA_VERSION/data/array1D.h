/***************************************************************************
 *
 * Authors:    Yuxuan Chen
 *
 * Copyright (C) 2021 Pattern Recognition and Bioinformatics Group, Shanghai Jiao Tong University
 *
 * Licensed under the MIT License (see LICENSE for details)
 *
 * All comments concerning this program package may be sent to the e-mail address 'yxchen11@sjtu.edu.cn'
 ***************************************************************************/

#ifndef ALIGNMENT_ARRAY1D_H
#define ALIGNMENT_ARRAY1D_H


#include <cmath>
#include <cstring>
#include <typeinfo>
#include <string>
#include <iostream>


class baseException : std::exception {
public:
    std::string msg;

public:
    explicit baseException(std::string s) : msg(std::move(s)) {
        std::cerr << msg << std::endl;
    }
};


///only support numeric type: int, double, float
template<typename T>
class arrayReal {
public:
    int maxLengthInit = 100;
    long long maxLength{};
    long long length{};
    T *data = nullptr;

public:
    arrayReal() {
        data = new T[maxLengthInit];
        length = 0;
        maxLength = maxLengthInit;
    }

    arrayReal(const arrayReal<T> &ary) {
        data = new T[ary.maxLength];
        length = ary.length;
        maxLength = ary.maxLength;
        memcpy(data, ary.data, length * sizeof(T));
    }

    explicit arrayReal(const arrayReal<T> *ary) : arrayReal<T>(*ary) {
    }

    arrayReal(const T *d, long long l) {
        maxLength = l > maxLengthInit ? l : maxLengthInit;
        data = new T[maxLength];
        length = l;
        memcpy(data, d, l * sizeof(T));
    }

    ///zeros
    explicit arrayReal(long long l) {
        maxLength = l > maxLengthInit ? l : maxLengthInit;
        data = new T[maxLength];
        for (long long i = 0; i < l; ++i)
            data[i] = 0;
        length = l;
    }

    ~arrayReal() {
        if (data) {
            delete[] data;
            data = nullptr;
        }
    }

    arrayReal(const std::initializer_list<T> &initList) {
        maxLength = initList.size() > maxLengthInit ? initList.size() : maxLengthInit;
        data = new T[maxLength];
        length = initList.size();
        const T *p = initList.begin();
        memcpy(data, p, length * sizeof(T));
    }

    void append(T value) {
        if (length + 1 > maxLength) {
            T *newData = new T[length + 1];
            memcpy(newData, data, length * sizeof(T));
            newData[length] = value;
            delete[] data;
            data = newData;
            length++;
            maxLength++;
        } else {
            data[length] = value;
            length++;
        }
    }

    void remove() {
        if (length > maxLengthInit) {
            maxLength = length - 1 > maxLengthInit ? length - 1 : maxLengthInit;
            T *newData = new T[maxLength];
            memcpy(newData, data, (length - 1) * sizeof(T));
            delete[] data;
            data = newData;
            length--;
        } else if (length == 0)
            throw baseException("Error: Array is empty!");
        else
            length--;
    }

    arrayReal<T> get(long long index, long long l) const {
        if (index + l > length)
            throw baseException("Error: Index out of range!");
        arrayReal<T> ary = arrayReal<T>(l);
        memcpy(ary.data, data + index, l * sizeof(T));
        return ary;
    }

    void set(long long index, const arrayReal<T> &ary) {
        if (index + ary.length > length)
            throw baseException("Error: Index out of range!");
        memcpy(data + index, ary.data, ary.length * sizeof(T));
    }

    arrayReal<int> isNotInf() const {
        arrayReal<int> index = arrayReal<int>(length);
        for (long long i = 0; i < length; ++i)
            index[i] = data[i] == INFINITY ? 0 : 1;
        return index;
    }

    T sum() const {
        T s = 0;
        for (long long i = 0; i < length; ++i) {
            s += data[i];
        }
        return s;
    }

    T max() const {
        T r = -INFINITY;
        for (long long i = 0; i < length; ++i)
            r = r < data[i] ? data[i] : r;
        return r;
    }

    long long argmax() const {
        T r = -INFINITY;
        long long index = 0;
        for (long long i = 0; i < length; ++i) {
            if (r < data[i]) {
                r = data[i];
                index = i;
            }
        }
        return index;
    }

    T min() const {
        T r = INFINITY;
        for (long long i = 0; i < length; ++i)
            r = r > data[i] ? data[i] : r;
        return r;
    }

    long long argmin() const {
        T r = INFINITY;
        long long index = 0;
        for (long long i = 0; i < length; ++i) {
            if (r > data[i]) {
                r = data[i];
                index = i;
            }
        }
        return index;
    }

    float mean() const {
        float s = 0;
        for (long long i = 0; i < length; ++i) {
            s += static_cast<double>(data[i]) / length;
        }
        return s;
    }

    float var() const {
        float m = mean();
        float var = 0;
        for (long long i = 0; i < length; ++i) {
            var += (data[i] - m) * (data[i] - m) / length;
        }
        return var;
    }

    arrayReal<T> copy() const {
        arrayReal<T> ary = arrayReal<T>(this);
        return ary;
    }

    bool operator==(const arrayReal<T> &ary) const {
        if (typeid(ary.data) != typeid(data) || length != ary.length)
            return false;
        for (long long i = 0; i < length; ++i)
            if (ary.data[i] != data[i])
                return false;
        return true;
    }

    bool equalList(const std::initializer_list<T> &initList) const {
        if (length != initList.size())
            return false;
        const T *p;
        long long i;
        for (p = initList.begin(), i = 0; p != initList.end(); p++, ++i)
            if (data[i] != *p)
                return false;
        return true;
    }

    T &operator[](long long index) const {
        return data[index];
    }

    arrayReal<T> &operator=(const arrayReal<T> &ary) {
        if (this == &ary) {
            return *this;
        }
        if (length == ary.length) {
            memcpy(data, ary.data, length * sizeof(T));
            return *this;
        } else {
            delete[] data;
            data = new T[ary.maxLength];
            length = ary.length;
            maxLength = ary.maxLength;
            memcpy(data, ary.data, length * sizeof(T));
            return *this;
        }
    }

    template<typename T1>
    arrayReal<T> operator+(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = ary[i] + value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator-(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = ary[i] - value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator*(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = ary[i] * value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator/(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = value == 0 ? INFINITY : ary[i] / value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator+(const arrayReal<T1> ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = aryOut[i] + ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator-(const arrayReal<T1> ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = aryOut[i] - ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator*(const arrayReal<T1> ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = aryOut[i] * ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator/(const arrayReal<T1> ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            auto r = ary[i] == 0 ? INFINITY : aryOut[i] / ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    bool empty() const {
        return length == 0;
    }

    arrayReal<T> abs() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i)
            ary[i] = std::abs(ary[i]);
        return ary;
    }

    arrayReal<T> log() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i)
            ary[i] = std::log(ary[i]);
        return ary;
    }

    arrayReal<T> sin() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i)
            ary[i] = std::sin(ary[i]);
        return ary;
    }

    arrayReal<T> cos() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i)
            ary[i] = std::cos(ary[i]);
        return ary;
    }

    arrayReal<T> sqrt() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i)
            ary[i] = std::sqrt(ary[i]);
        return ary;
    }

    template<typename T_IN>
    arrayReal<T_IN> asType() const {
        arrayReal<T_IN> ary = arrayReal<T_IN>(length);
        for (long long i = 0; i < length; ++i)
            ary[i] = static_cast<T_IN>(data[i]);
        return ary;
    }

    arrayReal<T> clip(T minValue, T maxValue) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < length; ++i) {
            if (ary.data[i] < minValue)
                ary.data[i] = minValue;
            if (ary.data[i] > maxValue)
                ary.data[i] = maxValue;
        }
        return ary;
    }

    arrayReal<T> sort() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 1; i < ary.length; ++i)
            for (long long j = i; j > 0; --j) {
                if (ary[j] < ary[j - 1]) {
                    T temp = ary[j];
                    ary[j] = ary[j - 1];
                    ary[j - 1] = temp;
                } else break;
            }
        return ary;
    }

    arrayReal<T> replace(T value, T valueR) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (long long i = 0; i < ary.length; ++i) {
            ary[i] = ary[i] == value ? valueR : ary[i];
        }
        return ary;
    }

    long long exist(T value) const {
        long long count = 0;
        for (long long i = 0; i < length; ++i) {
            if (data[i] == value)
                ++count;
        }
        return count;
    }

    arrayReal<T> unique() const {
        arrayReal<T> ary;
        for (long long i = 0; i < length; ++i) {
            if (ary.exist(data[i]) == 0)
                ary.append(data[i]);
        }
        return ary;
    }

    arrayReal<T> intersectionA(const arrayReal<T> &ary) const {
        arrayReal<T> temp = unique();
        arrayReal<T> aryT = ary.unique();
        arrayReal<T> r;
        for (long long i = 0; i < aryT.length; ++i) {
            if (temp.exist(aryT[i]) != 0)
                r.append(aryT[i]);
        }
        return r;
    }

    arrayReal<T> unionA(const arrayReal<T> &ary) const {
        arrayReal<T> r = unique();
        arrayReal<T> aryT = ary.unique();
        for (long long i = 0; i < aryT.length; ++i) {
            r.append(aryT[i]);
        }
        r = r.unique();
        return r;
    }

    void print() const {
        std::cout << typeid(T).name() << '\t';
        for (long long i = 0; i < length; ++i)
            std::cout << data[i] << ' ';
        std::cout << std::endl;
    }

};


/// real value in odd index and imaginary value in even index
class arrayComplex {
public:
    int maxLengthInit = 100;
    long long maxLength{};
    long long length{};
    double *data = nullptr;

public:
    arrayComplex();

    arrayComplex(const arrayComplex &ary);

    explicit arrayComplex(const arrayComplex *ary);

    arrayComplex(const double d[][2], long long l);

    ///zeros
    explicit arrayComplex(long long l);

    ~arrayComplex();

    void append(const double value[2]);

    void append(double value);

    void remove();

    double *get(long long index) const;

    arrayComplex get(long long index, long long l) const;

    void set(long long index, const arrayComplex &ary);

    void set(long long index, const double value[2]);

    arrayComplex copy() const;

    bool empty() const;

    arrayReal<double> abs() const;

    arrayReal<double> angle() const;

    arrayReal<double> real() const;

    arrayReal<double> imag() const;

    double *sum() const;

    double *mean() const;

    double *operator[](long long index) const;

    arrayComplex &operator=(const arrayComplex &ary);

    arrayComplex conj() const;

    void print() const;

    arrayComplex operator+(double value) const;

    arrayComplex operator-(double value) const;

    arrayComplex operator*(double value) const;

    arrayComplex operator/(double value) const;

    arrayComplex operator+(const double value[2]) const;

    arrayComplex operator-(const double value[2]) const;

    arrayComplex operator*(const double value[2]) const;

    arrayComplex operator/(const double value[2]) const;

    ///only judge on real value
    arrayReal<int> isNotInf() const;

    arrayComplex operator+(const arrayComplex &ary) const;

    arrayComplex operator-(const arrayComplex &ary) const;

    arrayComplex operator*(const arrayComplex &ary) const;

    arrayComplex operator/(const arrayComplex &ary) const;

    template<typename T>
    arrayComplex operator+(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (long long i = 0; i < length; ++i)
            aryOut.data[2 * i] += ary.data[i];
        return aryOut;
    }

    template<typename T>
    arrayComplex operator-(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (long long i = 0; i < length; ++i)
            aryOut.data[2 * i] -= ary.data[i];
        return aryOut;
    }

    template<typename T>
    arrayComplex operator*(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (long long i = 0; i < length; ++i) {
            aryOut.data[2 * i] *= ary.data[i];
            aryOut.data[2 * i + 1] *= ary.data[i];
        }
        return aryOut;
    }

    template<typename T>
    arrayComplex operator/(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (long long i = 0; i < length; ++i) {
            if (ary.data[i] == 0) {
                aryOut.data[2 * i] = INFINITY;
                aryOut.data[2 * i + 1] = 0;
            } else {
                aryOut.data[2 * i] /= ary.data[i];
                aryOut.data[2 * i + 1] /= ary.data[i];
            }
        }
        return aryOut;
    }

};


arrayComplex bind(const arrayReal<double> &abs, const arrayReal<double> &angle);


#endif //ALIGNMENT_ARRAY1D_H