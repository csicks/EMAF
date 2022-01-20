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

#ifndef EM_POLAR_H
#define EM_POLAR_H

/** @file
 * this file contains class used in XMIPP alignment which transform image to polar coordinates
*/

#include "image2D.h"
#include "fft.h"
#include <vector>

/** @brief
 * class used in XMIPP which transform image to polar coordinates
*/
class polar {
public:
    arrayReal<double> radius; // radius of polar rings
    std::vector<arrayReal<double>> rings;
    std::vector<arrayComplex> ringsFFT; // FFT of polar rings

public:
    polar();

    polar(const polar &p);

    explicit polar(const imageReal<double> &img);

    /// get rings between two radius
    void getPolar(const imageReal<double> &img, int r1, int r2);

    /// apply FFT to rings
    void fftRings();

    /// get conjugate of FFT rings
    polar conj() const;

    /// get mean value and variance value of rings
    double *statistics();

    /// normalize rings using formula v=(v-mean)/var
    void normalize(double mean, double var);
};

/// linear interpolation
inline double linear(double value, double left, double right);

/// bi-linear interpolation for image pixel
double biLinear(const imageReal<double> &img, double x, double y);

/// compute correlation of two images in polar coordinates
arrayReal<double> polarCorrelation(const imageReal<double> &img1, const imageReal<double> &img2);

/// get best angle from the maximum value of array (learned from correlation)
double bestAngle(const arrayReal<double> &ary);

#endif //EM_POLAR_H