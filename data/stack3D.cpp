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

#include "stack3D.h"

stackComplex::stackComplex()
{
    shape[0] = 0;
    shape[1] = 0;
    shape[2] = 0;
    data = arrayComplex();
}

stackComplex::stackComplex(const int *s)
{
    shape[0] = s[0];
    shape[1] = s[1];
    shape[2] = s[2];
    data = arrayComplex(shape[0] * shape[1] * shape[2]);
}

stackComplex::stackComplex(const stackComplex &stk)
{
    shape[0] = stk.shape[0];
    shape[1] = stk.shape[1];
    shape[2] = stk.shape[2];
    data = arrayComplex(stk.data);
}

double *stackComplex::get(int i, int j, int k) const
{
    return data[i * shape[1] * shape[2] + j * shape[2] + k];
}

void stackComplex::set(int i, int j, int k, const double value[2])
{
    data.set(i * shape[1] * shape[2] + j * shape[2] + k, value);
}

imageComplex stackComplex::pieceGet(int index) const
{
    int imgShape[2] = {shape[1], shape[2]};
    imageComplex img = imageComplex(imgShape);
    img.data.set(0, data.get(index * shape[1] * shape[2], shape[1] * shape[2]));
    return img;
}

void stackComplex::pieceSet(int index, const imageComplex &img)
{
    data.set(index * shape[1] * shape[2], img.data);
}

double *stackComplex::mean() const
{
    return data.mean();
}

bool stackComplex::empty() const
{
    return shape[0] == 0;
}

void stackComplex::append(const imageComplex &img)
{
    if (img.shape[0] != shape[1] || img.shape[1] != shape[2])
        throw baseException("Error: Input imageReal is not of the same shape of images in stackReal!");
    if (data.length + img.shape[0] * img.shape[1] > data.maxLength)
    {
        arrayComplex temp = data.copy();
        data = arrayComplex(data.length + img.shape[0] * img.shape[1]);
        memcpy(data.data, temp.data, temp.length * sizeof(double) * 2);
        memcpy(data.data + temp.length, img.data.data, img.shape[0] * img.shape[1] * sizeof(double) * 2);
    }
    else
    {
        memcpy(data.data + data.length, img.data.data, img.shape[0] * img.shape[1] * sizeof(double) * 2);
    }
    shape[0] += 1;
}

void stackComplex::remove()
{
    if (shape[0] == 0)
        throw baseException("Error: Stack is empty!");
    else if (data.length > data.maxLengthInit)
    {
        arrayComplex temp = data.copy();
        data.maxLength = data.length - shape[1] * shape[2] > data.maxLengthInit ? data.length - shape[1] * shape[2] : data.maxLengthInit;
        data = arrayComplex(data.maxLength);
        data.length = temp.length - shape[1] * shape[2];
        memcpy(data.data, temp.data, data.length * sizeof(double) * 2);
        shape[0] -= 1;
    }
    else
    {
        shape[0] -= 1;
        data.length -= shape[1] * shape[2];
    }
}

stackComplex &stackComplex::operator=(const stackComplex &stk)
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

stackComplex image2Stack(const imageComplex &img)
{
    int stkShape[3] = {1, img.shape[0], img.shape[1]};
    stackComplex stk = stackComplex(stkShape);
    stk.pieceSet(0, img);
    return stk;
}

stackComplex image2Stack(const imageComplex img[], int l)
{
    int stkShape[3] = {l, img[0].shape[0], img[0].shape[1]};
    for (int i = 0; i < l; ++i)
        if (img[i].shape[0] != stkShape[1] || img[i].shape[1] != stkShape[2])
            throw baseException("Error: Images in stackReal are not of the same shape!");
    stackComplex stk = stackComplex(stkShape);
    for (int i = 0; i < l; ++i)
        stk.pieceSet(i, img[i]);
    return stk;
}