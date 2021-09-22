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

#ifndef ALIGNMENT_POLAR_H
#define ALIGNMENT_POLAR_H

#include "image2D.h"
#include "fft.h"
#include <vector>

#define ACCURACY 1e-6

class polar
{
public:
    arrayReal<double> radius;
    std::vector<arrayReal<double>> rings;
    std::vector<arrayComplex> ringsFFT;

public:
    polar();

    polar(const polar &p);

    explicit polar(const imageReal<float> &img);

    void getPolar(const imageReal<float> &img, int r1, int r2);

    void fftRings();

    polar conj() const;

    double *statistics();

    void normalize(double mean, double var);
};

double linear(double value, double left, double right);

double biLinear(const imageReal<float> &img, double x, double y);

arrayReal<double> polarCorrelation(const imageReal<float> &img1, const imageReal<float> &img2);

double bestAngle(const arrayReal<double> &ary);

#endif //ALIGNMENT_POLAR_H