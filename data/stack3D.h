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

#ifndef EM_STACK3D_H
#define EM_STACK3D_H

/** @file
 * this file contains data structure for three-dimensional image stack (real and complex)
 * TODO: may combine real/complex cases using std::complex (recommend) or _Complex/__complex__ (not recommend)
 * TODO: may replace get/set with () and remove [] (which will be supported in C++23)
*/

#include "array1D.h"
#include "image2D.h"

/** @brief
 * class for three-dimensional real image stack, whose size is changeable
*/
template<typename T>
class stackReal {
public:
    int shape[3]{}; // shape of image stack
    arrayReal<T> data; // data stored in one-dimensional array

public:
    stackReal() {
        shape[0] = 0;
        shape[1] = 0;
        shape[2] = 0;
        data = arrayReal<T>();
    }

    explicit stackReal(const int s[3]) {
        shape[0] = s[0];
        shape[1] = s[1];
        shape[2] = s[2];
        data = arrayReal<T>(static_cast<size_t>(shape[0]) * shape[1] * shape[2]);
    }

    stackReal(const stackReal<T> &stk) {
        shape[0] = stk.shape[0];
        shape[1] = stk.shape[1];
        shape[2] = stk.shape[2];
        data = arrayReal<T>(stk.data);
    }

    /// @return value at position (i,j,k), i.e. S[i,j,k]
    T get(int i, int j, int k) const {
        return data[static_cast<size_t>(i) * shape[1] * shape[2] + static_cast<size_t>(j) * shape[2] + k];
    }

    /// set S[i,j,k] to value
    void set(int i, int j, int k, T value) {
        data[static_cast<size_t>(i) * shape[1] * shape[2] + static_cast<size_t>(j) * shape[2] + k] = value;
    }

    /// @return the ith piece of image stack
    imageReal<T> pieceGet(int index) const {
        int s[2] = {shape[1], shape[2]};
        imageReal<T> img = imageReal<T>(s);
        img.data.set(0, data.get(static_cast<size_t>(index) * shape[1] * shape[2],
                                 static_cast<size_t>(shape[1]) * shape[2]));
        return img;
    }

    /// set the ith piece of image stack to img
    void pieceSet(int index, const imageReal<T> &img) {
        data.set(static_cast<size_t>(index) * shape[1] * shape[2], img.data);
    }

    /// @return maximum value of image stack
    T max() const {
        return data.max();
    }

    /// @return minimum value of image stack
    T min() const {
        return data.min();
    }

    /// @return average value of image stack
    T mean() const {
        return data.mean();
    }

    /// @return whether image stack is empty
    bool empty() const {
        return shape[0] == 0;
    }

    /// append image to the end of image stack
    void append(const imageReal<T> &img) {
        if (img.shape[0] != shape[1] or img.shape[1] != shape[2])
            throw baseException("Error: Input imageReal is not of the same shape of images in stackReal!");
        if (data.length + static_cast<size_t>(img.shape[0]) * img.shape[1] > data.maxLength) {
            arrayReal<T> temp = data.copy();
            data = arrayReal<T>(data.length + static_cast<size_t>(img.shape[0]) * img.shape[1]);
            memcpy(data.data, temp.data, temp.length * sizeof(T));
            memcpy(data.data + temp.length, img.data.data,
                   static_cast<size_t>(img.shape[0]) * img.shape[1] * sizeof(T));
        } else {
            memcpy(data.data + data.length, img.data.data,
                   static_cast<size_t>(img.shape[0]) * img.shape[1] * sizeof(T));
        }
        shape[0] += 1;
    }

    /// remove image at the end of image stack
    void remove() {
        if (shape[0] == 0)
            throw baseException("Error: Stack is empty!");
        else if (data.length > data.maxLengthInit) {
            arrayReal<T> temp = data.copy();
            data.maxLength = data.length - static_cast<size_t>(shape[1]) * shape[2] > data.maxLengthInit ?
                             data.length - static_cast<size_t>(shape[1]) * shape[2] : data.maxLengthInit;
            data = arrayReal<T>(data.maxLength);
            data.length = temp.length - static_cast<size_t>(shape[1]) * shape[2];
            memcpy(data.data, temp.data, data.length * sizeof(T));
            shape[0] -= 1;
        } else {
            shape[0] -= 1;
            data.length -= static_cast<size_t>(shape[1]) * shape[2];
        }
    }

    stackReal<T> &operator=(const stackReal<T> &stk) {
        if (this == &stk) {
            return *this;
        }
        shape[0] = stk.shape[0];
        shape[1] = stk.shape[1];
        shape[2] = stk.shape[2];
        data = stk.data;
        return *this;
    }

    /// get pieces from index start to end, i.e. S[i:j,:,:]
    stackReal<T> pieces(int start, int end) {
        int s[3] = {end - start + 1, shape[1], shape[2]};
        stackReal<T> stk = stackReal<T>(s);
        for (int i = start; i < end + 1; ++i) {
            stk.pieceSet(i - start, pieceGet(i));
        }
        return stk;
    }

    /// @return average image of image stack
    imageReal<T> stackAvg() {
        int l = shape[0];
        int s[2] = {shape[1], shape[2]};
        imageReal<T> img = imageReal<T>(s);
        for (int i = 0; i < l; ++i) {
            img = img + pieceGet(i) / l;
        }
        return img;
    }

};

/// convert one image to image stack
template<typename T>
stackReal<T> image2Stack(const imageReal<T> &img) {
    int stkShape[3] = {1, img.shape[0], img.shape[1]};
    stackReal<T> stk = stackReal<T>(stkShape);
    stk.pieceSet(0, img);
    return stk;
}

/// convert one image array to image stack
template<typename T>
stackReal<T> image2Stack(const imageReal<T> img[], int l) {
    int stkShape[3] = {l, img[0].shape[0], img[0].shape[1]};
    for (int i = 0; i < l; ++i)
        if (img[i].shape[0] != stkShape[1] or img[i].shape[1] != stkShape[2])
            throw baseException("Error: Images in stackReal are not of the same shape!");
    stackReal<T> stk = stackReal<T>(stkShape);
    for (int i = 0; i < l; ++i)
        stk.pieceSet(i, img[i]);
    return stk;
}

/** @brief
 * class for three-dimensional complex image stack, whose size is changeable
 * functions/variables appeared before in imageReal are ignored
*/
class stackComplex {
public:
    int shape[3]{};
    arrayComplex data;

public:
    stackComplex();

    explicit stackComplex(const int s[3]);

    stackComplex(const stackComplex &stk);

    /// @attention remember to free memory
    double *get(int i, int j, int k) const;

    void set(int i, int j, int k, const double value[2]);

    imageComplex pieceGet(int index) const;

    void pieceSet(int index, const imageComplex &img);

    /// @attention remember to free memory
    double *mean() const;

    bool empty() const;

    void append(const imageComplex &img);

    void remove();

    stackComplex &operator=(const stackComplex &stk);

};

stackComplex image2Stack(const imageComplex &img);

stackComplex image2Stack(const imageComplex img[], int l);


#endif //EM_STACK3D_H
