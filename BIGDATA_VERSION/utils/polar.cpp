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

#include "polar.h"


polar::polar() {
    radius = arrayReal<double>();
    rings.clear();
    ringsFFT.clear();
}


polar::polar(const polar &p) {
    radius = p.radius;
    rings = p.rings;
    ringsFFT = p.ringsFFT;
}


polar::polar(const imageReal<float> &img) {
    radius = arrayReal<double>();
    rings.clear();
    ringsFFT.clear();
    int rS = img.shape[0] / 5;
    int rE = img.shape[0] / 2;
    getPolar(img, rS, rE);
    double *s = statistics();
    normalize(s[0], s[1]);
    delete[] s;
    fftRings();
}


void polar::getPolar(const imageReal<float> &img, int r1, int r2) {
    float centerX = img.shape[0] / 2;
    float centerY = img.shape[1] / 2;
    radius = arrayReal<double>(r2 - r1 + 1);

    for (int r = r1; r < r2 + 1; ++r) {
        int circle = static_cast<int>(2 * M_PI * 0.5 * r) * 2;
        double dAng = 2 * M_PI / circle;
        arrayReal<double> ring = arrayReal<double>(circle);
        for (int i = 0; i < circle; ++i) {
            double ang = i * dAng;
            double x = -std::sin(ang) * r + centerX;
            double y = std::cos(ang) * r + centerY;
            if (-1 <= x && x < 0)
                x = 0;
            if (img.shape[0] - 1 < x && x <= img.shape[0])
                x = img.shape[0] - 1;
            if (-1 <= y && y < 0)
                y = 0;
            if (img.shape[1] - 1 < y && y <= img.shape[1])
                y = img.shape[1] - 1;
            ring[i] = biLinear(img, x, y);
        }
        rings.push_back(ring);
        radius[r - r1] = r;
    }
}


void polar::fftRings() {
    for (int r = radius[0]; r < radius[radius.length - 1]; ++r) {
        arrayComplex ary = fftHalf(rings[r - radius[0]]);
        ringsFFT.push_back(ary);
    }
}


polar polar::conj() const {
    polar p = *this;
    for (auto &i : p.ringsFFT)
        i = i.conj();
    return p;
}


double *polar::statistics() {
    auto *s = new double[2];
    double sum = 0, sum2 = 0, N = 0;
    for (int i = 0; i < radius.length; ++i) {
        double r = radius[i];
        arrayReal<double> ring = rings[i];
        double weight = 2 * M_PI * r / ring.length;
        for (int j = 0; j < ring.length; ++j) {
            double value = ring[j];
            sum += value * weight;
            sum2 += value * value * weight;
            N += weight;
        }
    }
    sum2 /= N;
    double mean = sum / N;
    double stddev = std::sqrt(std::fabs(sum2 - mean * mean));
    s[0] = mean;
    s[1] = stddev;
    return s;
}


void polar::normalize(double mean, double var) {
    for (int i = 0; i < radius.length; ++i) {
        arrayReal<double> &ring = rings[i];
        for (int j = 0; j < ring.length; ++j)
            ring[j] = (ring[j] - mean) / var;
    }
}


double linear(double value, double left, double right) {
    return left + (right - left) * value;
}


double biLinear(const imageReal<float> &img, double x, double y) {
    int xl = std::floor(x);
    int xr = xl + 1;
    int yl = std::floor(y);
    int yr = yl + 1;
    double fx = x - xl;
    double fy = y - yl;
    double neighbour[4];
    if (0 <= xl && xl < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour[0] = img.get(xl, yl);
    else
        neighbour[0] = 0;
    if (0 <= xr && xr < img.shape[0] && 0 <= yl && yl < img.shape[1])
        neighbour[1] = img.get(xr, yl);
    else
        neighbour[1] = 0;
    if (0 <= xl && xl < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour[2] = img.get(xl, yr);
    else
        neighbour[2] = 0;
    if (0 <= xr && xr < img.shape[0] && 0 <= yr && yr < img.shape[1])
        neighbour[3] = img.get(xr, yr);
    else
        neighbour[3] = 0;
    double temp0 = linear(fx, neighbour[0], neighbour[1]);
    double temp1 = linear(fx, neighbour[2], neighbour[3]);
    double r = linear(fy, temp0, temp1);
    return r;
}


arrayReal<double> polarCorrelation(const imageReal<float> &img1, const imageReal<float> &img2) {
    if (img1.shape[0] != img2.shape[0] || img1.shape[1] != img2.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    polar p1 = polar(img1);
    polar p2 = polar(img2);
    int length = p1.ringsFFT.size();
    int maxSize = p1.ringsFFT[length - 1].length;
    arrayComplex temp = arrayComplex(maxSize);
    for (int i = 0; i < length; ++i) {
        arrayComplex dot = p1.ringsFFT[i] * p2.ringsFFT[i].conj();
        double circle = 2 * M_PI * p1.radius[i];
        dot = dot * circle;
        int maxIndex = p1.ringsFFT[i].length;
        for (int j = 0; j < maxIndex; ++j) {
            temp.data[2 * j] += dot.data[2 * j];
            temp.data[2 * j + 1] += dot.data[2 * j + 1];
        }
    }
    arrayReal<double> iary = ifftHalf(temp);
    return iary;
}


double bestAngle(const arrayReal<double> &ary) {
    int maxIndex = ary.argmax();
    double ang = static_cast<double>(maxIndex) / ary.length * 360.0;
    return ang;
}