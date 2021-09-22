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

#ifndef ALIGNMENT_FFT_H
#define ALIGNMENT_FFT_H

/** reference website:
           * http://www.fftw.org/fftw3_doc/ */

#include <fftw3.h>
#include <cstring>
#include <iostream>
#include "array1D.h"
#include "image2D.h"

template <typename T>
arrayComplex fft(const arrayReal<T> &ary)
{
    return fft(ary.template asType<double>());
}

arrayComplex fft(const arrayReal<double> &ary);

arrayComplex fft(const arrayComplex &ary);

arrayComplex ifftC(const arrayComplex &ary);

arrayReal<double> ifftR(const arrayComplex &ary);

arrayComplex fftShift(const arrayComplex &ary);

arrayComplex ifftShift(const arrayComplex &ary);

template <typename T>
arrayComplex fftHalf(const arrayReal<T> &ary)
{
    return fftHalf(ary.template asType<double>());
}

arrayComplex fftHalf(const arrayReal<double> &ary);

arrayReal<double> ifftHalf(const arrayComplex &ary);

template <typename T>
imageComplex fft2(const imageReal<T> &img)
{
    return fft2(img.template asType<double>());
}

imageComplex fft2(const imageReal<double> &img);

imageComplex fft2(const imageComplex &img);

imageComplex ifft2C(const imageComplex &img);

imageReal<double> ifft2R(const imageComplex &img);

imageComplex fftShift2(const imageComplex &img);

imageComplex ifftShift2(const imageComplex &img);

template <typename T>
imageComplex fft2Half(const imageReal<T> &img)
{
    return fft2Half(img.template asType<double>());
}

imageComplex fft2Half(const imageReal<double> &img);

imageReal<double> ifft2Half(const imageComplex &img);

#endif //ALIGNMENT_FFT_H