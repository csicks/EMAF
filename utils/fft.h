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

#ifndef EM_FFT_H
#define EM_FFT_H

/** @file
 * this file contains functions for performing one-dimensional, two-dimensional and three-dimensional FFTs
 * reference website: http://www.fftw.org/fftw3_doc
 * TODO: improve the fftshift functions
 */

#include <fftw3.h>
#include <cstring>
#include <iostream>
#include "array1D.h"
#include "image2D.h"

/// one-dimensional FFT for real array
template <typename T>
arrayComplex fft(const arrayReal<T> &ary)
{
    return fft(ary.template asType<double>());
}

/// one-dimensional FFT for real array
arrayComplex fft(const arrayReal<double> &ary);

/// one-dimensional FFT for complex array
arrayComplex fft(const arrayComplex &ary);

/** one-dimensional inverse FFT
 * @return complex array
 */
arrayComplex ifftC(const arrayComplex &ary);

/** one-dimensional inverse FFT
 * @return real array
 */
arrayReal<double> ifftR(const arrayComplex &ary);

/// shift low frequencies to the center of array
arrayComplex fftShift(const arrayComplex &ary);

/// shift low frequencies back to the corner of array
arrayComplex ifftShift(const arrayComplex &ary);

/// one-dimensional FFT for real array, with only half the memory/space used
template <typename T>
arrayComplex fftHalf(const arrayReal<T> &ary)
{
    return fftHalf(ary.template asType<double>());
}

/// one-dimensional FFT for real array, with only half the memory/space used
arrayComplex fftHalf(const arrayReal<double> &ary);

/// one-dimensional inverse FFT for FFT only use half the memory/space
arrayReal<double> ifftHalf(const arrayComplex &ary, int originLength);

/// two-dimensional FFT for real image
template <typename T>
imageComplex fft2(const imageReal<T> &img)
{
    return fft2(img.template asType<double>());
}

/// two-dimensional FFT for real image
imageComplex fft2(const imageReal<double> &img);

/// two-dimensional FFT for complex image
imageComplex fft2(const imageComplex &img);

/** two-dimensional inverse FFT
 * @return complex image
 */
imageComplex ifft2C(const imageComplex &img);

/** two-dimensional inverse FFT
 * @return real image
 */
imageReal<double> ifft2R(const imageComplex &img);

/// shift low frequencies to the center of image
imageComplex fftShift2(const imageComplex &img);

/// shift low frequencies back to the corner of image
imageComplex ifftShift2(const imageComplex &img);

/// two-dimensional FFT for real image, with only half the memory/space used
template <typename T>
imageComplex fft2Half(const imageReal<T> &img)
{
    return fft2Half(img.template asType<double>());
}

/// two-dimensional FFT for real image, with only half the memory/space used
imageComplex fft2Half(const imageReal<double> &img);

/// two-dimensional inverse FFT for FFT only use half the memory/space
imageReal<double> ifft2Half(const imageComplex &img, int originLength);

#endif // EM_FFT_H