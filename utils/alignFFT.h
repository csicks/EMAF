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

#ifndef ALIGNMENT_ALIGNFFT_H
#define ALIGNMENT_ALIGNFFT_H


#include "image2D.h"
#include "fft.h"
#include "stack3D.h"
#include "polar.h"
#include "fileStack.h"


#define ALIGN_TH 1


struct phi {
    double ang = 0;
    double x = 0;
    double y = 0;

    phi operator+(phi p) const {
        phi r;
        r.ang = ang + p.ang;
        r.x = x + p.x;
        r.y = y + p.y;
        return r;
    }

    phi() = default;

    phi(const std::initializer_list<double> &initList) {
        const double *p = initList.begin();
        ang = *p++;
        x = *p++;
        y = *p++;
    }

};

phi shiftFFT(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi rotateFFT(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi shiftFFTX(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi rotateFFTX(const imageReal<float> &imgX, const imageReal<float> &imgC);

imageComplex polarFFT(const imageReal<float> &img);

phi alignFFT(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi angAll(double ang0, const imageReal<float> &imgX, const imageReal<float> &imgC, bool fineFlag = false);

phi alignFFT2(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi alignShiftPolar(const imageReal<float> &imgX, const imageReal<float> &imgC, int maxShiftX = 4, int maxShiftY = 4);

phi alignFFTX(const imageReal<float> &imgX, const imageReal<float> &imgC);

imageReal<float> applyPhi(const imageReal<float> &img, phi p);

imageReal<float> inversePhi(const imageReal<float> &img, phi p);

imageReal<float> reference(const stackReal<float> &stk, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC));

imageReal<float> reference(const stackReal<float> &stk, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width, bool fineFlag), int width,
                           bool fineFlag = false);

imageReal<float> reference(const stackReal<float> &stk, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int outWidth, int inWidth, float th, bool fineFlag),
                           int outWidth, int inWidth = -1, float th = -1, bool fineFlag = false);

imageReal<float> reference(const fileStack &stk, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC));

imageReal<float> centerImage(const imageReal<float> &img);

imageReal<float> initReference(const stackReal<float> &stk);

imageReal<float> initReference(const fileStack &stk);

void refinement(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), bool updateFlag=true);

void refinement(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width, bool fineFlag), int width,
                bool fineFlag = false, bool updateFlag=true);

void refinement(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int outWidth, int inWidth, float th,
                                                                       bool fineFlag), int outWidth, int inWidth = -1, float th = -1, bool fineFlag = false, bool updateFlag=true);

void outInfo(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC),
             std::vector<std::string> info, const std::string &xmdName);

void outInfo(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width, bool fineFlag),
             std::vector<std::string> info, const std::string &xmdName, int width, bool fineFlag = false);

void outInfo(stackReal<float> &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int outWidth, int inWidth, float th, bool fineFlag),
             std::vector<std::string> info, const std::string &xmdName, int outWidth, int inWidth = -1, float th = -1, bool fineFlag = false);

void refinement(fileStack &stk, imageReal<float> &ref, phi func(const imageReal<float> &imgX, const imageReal<float> &imgC));

void refinementWithAux(stackReal<float> &stk, stackReal<float> &stkAux, imageReal<float> &ref);

phi alignPower(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi alignPeak(const imageReal<float> &imgX, const imageReal<float> &imgC);

int circleEstimatePeak(const stackReal<float> &stk, bool fineFlag = false);

int circleEstimateShape(const stackReal<float> &stk, int inWidth = -1, float th = -1, bool fineFlag = false);

phi alignPeak(const imageReal<float> &imgX, const imageReal<float> &imgC, int width, bool fineFlag = false);

phi alignShape(const imageReal<float> &imgX, const imageReal<float> &imgC);

phi alignShape(const imageReal<float> &imgX, const imageReal<float> &imgC, int outWidth, int inWidth = -1, float th = -1, bool fineFlag = false);

phi alignParticle(const imageReal<float> &imgX, const imageReal<float> &imgC);

#endif //ALIGNMENT_ALIGNFFT_H