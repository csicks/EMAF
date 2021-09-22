/***************************************************************************
 *
 * Authors:    Yuxuan Chen
 *
 * Copyright (C) 2021 Pattern Recognition and Bioinformatics Group, Shanghai Jiao Tong University
 *
 * Licensed under the GNU General Public License v3.0 (see LICENSE for details)
 *
 * All comments concerning this program package may be sent to the e-mail address 'yxchen11@sjtu.edu.cn'
 ***************************************************************************/

#ifndef ALIGNMENT_STACK3D_H
#define ALIGNMENT_STACK3D_H

#include "array1D.h"
#include "image2D.h"

template <typename T>
class stackReal
{
public:
    int shape[3]{};
    arrayReal<T> data;

public:
    stackReal()
    {
        shape[0] = 0;
        shape[1] = 0;
        shape[2] = 0;
        data = arrayReal<T>();
    }

    explicit stackReal(const int *s)
    {
        shape[0] = s[0];
        shape[1] = s[1];
        shape[2] = s[2];
        data = arrayReal<T>(shape[0] * shape[1] * shape[2]);
    }

    stackReal(const stackReal<T> &stk)
    {
        shape[0] = stk.shape[0];
        shape[1] = stk.shape[1];
        shape[2] = stk.shape[2];
        data = arrayReal<T>(stk.data);
    }

    T get(int i, int j, int k) const
    {
        return data[i * shape[1] * shape[2] + j * shape[2] + k];
    }

    void set(int i, int j, int k, T value)
    {
        data[i * shape[1] * shape[2] + j * shape[2] + k] = value;
    }

    imageReal<T> pieceGet(int index) const
    {
        int imgShape[2] = {shape[1], shape[2]};
        imageReal<T> img = imageReal<T>(imgShape);
        img.data.set(0, data.get(index * shape[1] * shape[2], shape[1] * shape[2]));
        return img;
    }

    void pieceSet(int index, const imageReal<T> &img)
    {
        data.set(index * shape[1] * shape[2], img.data);
    }

    T max() const
    {
        return data.max();
    }

    T min() const
    {
        return data.min();
    }

    T mean() const
    {
        return data.mean();
    }

    bool empty() const
    {
        return shape[0] == 0;
    }

    void append(const imageReal<T> &img)
    {
        if (img.shape[0] != shape[1] || img.shape[1] != shape[2])
            throw baseException("Error: Input imageReal is not of the same shape of images in stackReal!");
        if (data.length + img.shape[0] * img.shape[1] > data.maxLength)
        {
            arrayReal<T> temp = data.copy();
            data = arrayReal<T>(data.length + img.shape[0] * img.shape[1]);
            memcpy(data.data, temp.data, temp.length * sizeof(T));
            memcpy(data.data + temp.length, img.data.data, img.shape[0] * img.shape[1] * sizeof(T));
        }
        else
        {
            memcpy(data.data + data.length, img.data.data, img.shape[0] * img.shape[1] * sizeof(T));
        }
        shape[0] += 1;
    }

    void remove()
    {
        if (shape[0] == 0)
            throw baseException("Error: Stack is empty!");
        else if (data.length > data.maxLengthInit)
        {
            arrayReal<T> temp = data.copy();
            data.maxLength = data.length - shape[1] * shape[2] > data.maxLengthInit ? data.length - shape[1] * shape[2] : data.maxLengthInit;
            data = arrayReal<T>(data.maxLength);
            data.length = temp.length - shape[1] * shape[2];
            memcpy(data.data, temp.data, data.length * sizeof(T));
            shape[0] -= 1;
        }
        else
        {
            shape[0] -= 1;
            data.length -= shape[1] * shape[2];
        }
    }

    stackReal<T> &operator=(const stackReal<T> &stk)
    {
        if (this == &stk)
        {
            return *this;
        }
        shape[0] = stk.shape[0];
        shape[1] = stk.shape[1];
        shape[2] = stk.shape[2];
        data = stk.data;
        return *this;
    }

    stackReal<T> pieces(int start, int end)
    {
        int s[3] = {end - start + 1, shape[1], shape[2]};
        stackReal<T> stk = stackReal<T>(s);
        for (int i = start; i < end + 1; ++i)
        {
            stk.pieceSet(i - start, pieceGet(i));
        }
        return stk;
    }

    imageReal<T> stackAvg()
    {
        int l = shape[0];
        int s[2] = {shape[1], shape[2]};
        imageReal<T> img = imageReal<T>(s);
        for (int i = 0; i < l; ++i)
        {
            img = img + pieceGet(i) / l;
        }
        return img;
    }
};

template <typename T>
stackReal<T> image2Stack(const imageReal<T> &img)
{
    int stkShape[3] = {1, img.shape[0], img.shape[1]};
    stackReal<T> stk = stackReal<T>(stkShape);
    stk.pieceSet(0, img);
    return stk;
}

template <typename T>
stackReal<T> image2Stack(const imageReal<T> img[], int l)
{
    int stkShape[3] = {l, img[0].shape[0], img[0].shape[1]};
    for (int i = 0; i < l; ++i)
        if (img[i].shape[0] != stkShape[1] || img[i].shape[1] != stkShape[2])
            throw baseException("Error: Images in stackReal are not of the same shape!");
    stackReal<T> stk = stackReal<T>(stkShape);
    for (int i = 0; i < l; ++i)
        stk.pieceSet(i, img[i]);
    return stk;
}

class stackComplex
{
public:
    int shape[3]{};
    arrayComplex data;

public:
    stackComplex();

    explicit stackComplex(const int *s);

    stackComplex(const stackComplex &stk);

    double *get(int i, int j, int k) const;

    void set(int i, int j, int k, const double value[2]);

    imageComplex pieceGet(int index) const;

    void pieceSet(int index, const imageComplex &img);

    double *mean() const;

    bool empty() const;

    void append(const imageComplex &img);

    void remove();

    stackComplex &operator=(const stackComplex &stk);
};

stackComplex image2Stack(const imageComplex &img);

stackComplex image2Stack(const imageComplex img[], int l);

#endif //ALIGNMENT_STACK3D_H