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

#ifndef EM_ARRAY1D_H
#define EM_ARRAY1D_H

/** @file
 * this file contains data structure for one-dimensional numeric array (real and complex)
 * TODO: may combine real/complex cases using std::complex (recommend) or _Complex/__complex__ (not recommend)
 * @code
 * complex<double> f[100] = {c, c, c,};
 * auto out_f = reinterpret_cast<double (&)[][2] > (f);
 * double tf[100][2];
 * memcpy(tf, out_f, 200 * sizeof(double));
 * @endcode
 * TODO: may replace [] with () which would unify the style of all dimensional data
*/

#include <cmath>
#include <cstring>
#include <typeinfo>
#include <iostream>
#include "logger.h"

/** @brief
 * class for one-dimensional array, which only support numeric type like int, float and double
 * the size is changeable but additional cost would be required
*/
template<typename T>
class arrayReal {
public:
    size_t maxLengthInit = 100; // initial maximum size of array, used for small arrays
    size_t maxLength{}; // current maximum size of array
    size_t length{}; // current length of array
    T *data = nullptr; // pointer to data

public:
    arrayReal() {
        length = 0;
        maxLength = maxLengthInit;
        data = new T[maxLengthInit];
    }

    arrayReal(const arrayReal<T> &ary) {
        length = ary.length;
        maxLength = ary.maxLength;
        data = new T[ary.maxLength];
        memcpy(data, ary.data, length * sizeof(T));
    }

    explicit arrayReal(const arrayReal<T> *ary) : arrayReal<T>(*ary) {
    }

    arrayReal(const T *d, size_t l) {
        length = l;
        maxLength = l > maxLengthInit ? l : maxLengthInit;
        data = new T[maxLength];
        memcpy(data, d, l * sizeof(T));
    }

    /// initialization of length l and filled with zero values
    explicit arrayReal(size_t l) {
        length = l;
        maxLength = l > maxLengthInit ? l : maxLengthInit;
        data = new T[maxLength];
        for (size_t i = 0; i < l; ++i)
            data[i] = 0;
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

    /// append value at the end of array
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

    /// remove value at the end of array
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

    /// get values from index to index+l
    arrayReal<T> get(size_t index, size_t l) const {
        if (index + l > length)
            throw baseException("Error: Index out of range!");
        arrayReal<T> ary = arrayReal<T>(l);
        memcpy(ary.data, data + index, l * sizeof(T));
        return ary;
    }

    /// set values to ary from index
    void set(size_t index, const arrayReal<T> &ary) {
        if (index + ary.length > length)
            throw baseException("Error: Index out of range!");
        memcpy(data + index, ary.data, ary.length * sizeof(T));
    }

    /** judge every value of array to be inf or not
     * 0 for not inf and 1 for inf
     * @return array filled with value 0 or 1
    */
    arrayReal<int> isNotInf() const {
        arrayReal<int> index = arrayReal<int>(length);
        for (size_t i = 0; i < length; ++i)
            index[i] = data[i] == INFINITY ? 0 : 1;
        return index;
    }

    /// @return sum of array
    T sum() const {
        T s = 0;
        for (size_t i = 0; i < length; ++i) {
            s += data[i];
        }
        return s;
    }

    /// @return maximum value of array
    T max() const {
        T r = -INFINITY;
        for (size_t i = 0; i < length; ++i)
            r = r < data[i] ? data[i] : r;
        return r;
    }

    /// @return index of maximum value in array
    size_t argmax() const {
        T r = -INFINITY;
        size_t index = 0;
        for (size_t i = 0; i < length; ++i) {
            if (r < data[i]) {
                r = data[i];
                index = i;
            }
        }
        return index;
    }

    /// @return minimum value of array
    T min() const {
        T r = INFINITY;
        for (size_t i = 0; i < length; ++i)
            r = r > data[i] ? data[i] : r;
        return r;
    }

    /// @return index of minimum value in array
    size_t argmin() const {
        T r = INFINITY;
        size_t index = 0;
        for (size_t i = 0; i < length; ++i) {
            if (r > data[i]) {
                r = data[i];
                index = i;
            }
        }
        return index;
    }

    /// @return average value of array
    double mean() const {
        double s = 0;
        for (size_t i = 0; i < length; ++i) {
            s += static_cast<double>(data[i]) / length;
        }
        return s;
    }

    /// @return variance of array
    double var() const {
        double m = mean();
        double var = 0;
        for (size_t i = 0; i < length; ++i) {
            var += (data[i] - m) * (data[i] - m) / length;
        }
        return var;
    }

    /// @return same array of this one which doesn't share memory
    arrayReal<T> copy() const {
        arrayReal<T> ary = arrayReal<T>(this);
        return ary;
    }

    bool operator==(const arrayReal<T> &ary) const {
        if (typeid(ary.data) != typeid(data) or length != ary.length)
            return false;
        for (size_t i = 0; i < length; ++i)
            if (ary.data[i] != data[i])
                return false;
        return true;
    }

    /// @return whether array have the same value with the input list
    bool equalList(const std::initializer_list<T> &initList) const {
        if (length != initList.size())
            return false;
        const T *p;
        size_t i;
        for (p = initList.begin(), i = 0; p != initList.end(); p++, ++i)
            if (data[i] != *p)
                return false;
        return true;
    }

    T &operator[](size_t index) const {
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
        for (size_t i = 0; i < length; ++i) {
            auto r = ary[i] + value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator-(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = ary[i] - value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator*(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = ary[i] * value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator/(T1 value) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = value == 0 ? INFINITY : ary[i] / value;
            ary[i] = static_cast<T>(r);
        }
        return ary;
    }

    template<typename T1>
    arrayReal<T> operator+(const arrayReal<T1> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = aryOut[i] + ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator-(const arrayReal<T1> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = aryOut[i] - ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator*(const arrayReal<T1> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = aryOut[i] * ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    template<typename T1>
    arrayReal<T> operator/(const arrayReal<T1> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayReal<T> aryOut = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            auto r = ary[i] == 0 ? INFINITY : aryOut[i] / ary[i];
            aryOut[i] = static_cast<T>(r);
        }
        return aryOut;
    }

    /// @return whether array is emtpy
    bool empty() const {
        return length == 0;
    }

    /// @return absolute value of array
    arrayReal<T> abs() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i)
            ary[i] = std::abs(ary[i]);
        return ary;
    }

    /// @return logarithm_e value of array
    arrayReal<T> log() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i)
            ary[i] = std::log(ary[i]);
        return ary;
    }

    /// @return sine value of array
    arrayReal<T> sin() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i)
            ary[i] = std::sin(ary[i]);
        return ary;
    }

    /// @return cosine value of array
    arrayReal<T> cos() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i)
            ary[i] = std::cos(ary[i]);
        return ary;
    }

    /// @return square root value of array
    arrayReal<T> sqrt() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i)
            ary[i] = std::sqrt(ary[i]);
        return ary;
    }

    /** convert data type to another one
     * @example ary.asType<double>()
    */
    template<typename T1>
    arrayReal<T1> asType() const {
        arrayReal<T1> ary = arrayReal<T1>(length);
        for (size_t i = 0; i < length; ++i)
            ary[i] = static_cast<T1>(data[i]);
        return ary;
    }

    /** force the value range of array to be (minValue,maxValue)
     * value larger than maxValue will be set to maxValue
     * value smaller than minValue will be set to maxValue
    */
    arrayReal<T> clip(T minValue, T maxValue) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < length; ++i) {
            if (ary.data[i] < minValue)
                ary.data[i] = minValue;
            if (ary.data[i] > maxValue)
                ary.data[i] = maxValue;
        }
        return ary;
    }

    /// sort array in increasing order
    arrayReal<T> sort() const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 1; i < ary.length; ++i)
            for (size_t j = i; j > 0; --j) {
                if (ary[j] < ary[j - 1]) {
                    T temp = ary[j];
                    ary[j] = ary[j - 1];
                    ary[j - 1] = temp;
                } else
                    break;
            }
        return ary;
    }

    /// replace every value in array with valueR
    arrayReal<T> replace(T value, T valueR) const {
        arrayReal<T> ary = arrayReal<T>(this);
        for (size_t i = 0; i < ary.length; ++i) {
            ary[i] = ary[i] == value ? valueR : ary[i];
        }
        return ary;
    }

    /// @return number of value in array
    size_t exist(T value) const {
        size_t count = 0;
        for (size_t i = 0; i < length; ++i) {
            if (data[i] == value)
                ++count;
        }
        return count;
    }

    /// @return unique value in array
    arrayReal<T> unique() const {
        arrayReal<T> ary;
        for (size_t i = 0; i < length; ++i) {
            if (ary.exist(data[i]) == 0)
                ary.append(data[i]);
        }
        return ary;
    }

    /// @return intersection of two arrays
    arrayReal<T> intersectionA(const arrayReal<T> &ary) const {
        arrayReal<T> temp = unique();
        arrayReal<T> aryT = ary.unique();
        arrayReal<T> r;
        for (size_t i = 0; i < aryT.length; ++i) {
            if (temp.exist(aryT[i]) != 0)
                r.append(aryT[i]);
        }
        return r;
    }

    /// @return union of two arrays
    arrayReal<T> unionA(const arrayReal<T> &ary) const {
        arrayReal<T> r = unique();
        arrayReal<T> aryT = ary.unique();
        for (size_t i = 0; i < aryT.length; ++i) {
            r.append(aryT[i]);
        }
        r = r.unique();
        return r;
    }

    /// print array, with data type first
    void print() const {
        std::cout << typeid(T).name() << '\t';
        for (size_t i = 0; i < length; ++i)
            std::cout << data[i] << ' ';
        std::cout << std::endl;
    }
};

/** @brief
 * class for one-dimensional array, which support complex value
 * the complex class in c++ is not used
 * instead, real value is stored in odd index and imaginary value is stored in even index
 * functions/variables appeared before in arrayReal are ignored
*/
class arrayComplex {
public:
    size_t maxLengthInit = 100;
    size_t maxLength{};
    size_t length{};
    double *data = nullptr;

public:
    arrayComplex();

    arrayComplex(const arrayComplex &ary);

    explicit arrayComplex(const arrayComplex *ary);

    arrayComplex(const double d[][2], size_t l);

    explicit arrayComplex(size_t l);

    ~arrayComplex();

    void append(const double value[2]);

    void append(double value);

    void remove();

    double *get(size_t index) const;

    arrayComplex get(size_t index, size_t l) const;

    void set(size_t index, const arrayComplex &ary);

    void set(size_t index, const double value[2]);

    arrayComplex copy() const;

    bool empty() const;

    arrayReal<double> abs() const;

    /// @return complex angle of complex array
    arrayReal<double> angle() const;

    /// @return real part of complex array
    arrayReal<double> real() const;

    /// @return imaginary part of complex array
    arrayReal<double> imag() const;

    /// @attention remember to free memory
    double *sum() const;

    /// @attention remember to free memory
    double *mean() const;

    /** @attention remember to free memory
     * @attention assignment value, i.e. array[i]=v, doesn't work!
    */
    double *operator[](size_t index) const;

    arrayComplex &operator=(const arrayComplex &ary);

    /// @return conjugate of complex array
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

    /// only judge on real value
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
        for (size_t i = 0; i < length; ++i)
            aryOut.data[2 * i] += ary.data[i];
        return aryOut;
    }

    template<typename T>
    arrayComplex operator-(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (size_t i = 0; i < length; ++i)
            aryOut.data[2 * i] -= ary.data[i];
        return aryOut;
    }

    template<typename T>
    arrayComplex operator*(const arrayReal<T> &ary) const {
        if (ary.length != length)
            throw baseException("Error: Two arrays are not of the same shape!");
        arrayComplex aryOut = arrayComplex(this);
        for (size_t i = 0; i < length; ++i) {
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
        for (size_t i = 0; i < length; ++i) {
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

/// recover complex array according to its absolute value and complex angle
arrayComplex bind(const arrayReal<double> &abs, const arrayReal<double> &angle);

#endif //EM_ARRAY1D_H