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

#include "align.h"
#include "mrc.h"
#include "noise.h"
#include "particle.h"
#include "xmippXMD.h"

phi shiftFFT(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    if (imgX.shape[0] != imgC.shape[0] or imgX.shape[1] != imgC.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    imageComplex fx = fft2(imgX);
    imageComplex fc = fft2(imgC);
    imageComplex phaseF = fx * fc.conj() / (fx * fc.conj()).abs();
    imageReal<double> phase = ifft2R(phaseF);
    int *pos = phase.argmax();
    if (pos[0] > imgC.shape[0] / 2)
        pos[0] -= imgC.shape[0];
    if (pos[1] > imgC.shape[1] / 2)
        pos[1] -= imgC.shape[1];
    phi r = {0, static_cast<double>(pos[0]), static_cast<double>(pos[1])};
    delete[] pos;
    return r;
}

phi rotateFFT(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    if (imgX.shape[0] != imgC.shape[0] or imgX.shape[1] != imgC.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    imageReal<double> px = polarBack(imgX, false);
    imageReal<double> pc = polarBack(imgC, false);
    imageComplex fx = fft2(px);
    imageComplex fc = fft2(pc);
    imageComplex phaseF = fx * fc.conj() / (fx * fc.conj()).abs();
    imageReal<double> phase = ifft2R(phaseF);
    int *pos = phase.argmax();
    double angle = static_cast<double>(pos[0]) / imgC.shape[0] * 360;
    delete[] pos;
    phi r = {angle, 0, 0};
    return r;
}

phi shiftFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    if (imgX.shape[0] != imgC.shape[0] or imgX.shape[1] != imgC.shape[1])
        throw baseException("Error: Two images are not of the same shape!");
    int N = imgX.shape[0] * imgX.shape[1];

    imageComplex fx = fft2Half(imgX);
    imageComplex fc = fft2Half(imgC);
    imageComplex phaseF = fc * fx.conj() * N;
    imageReal<double> phase = ifft2Half(phaseF, imgX.shape[1]);
    int *pos = phase.argmax();
    if (pos[0] > imgC.shape[0] / 2)
        pos[0] -= imgC.shape[0];
    if (pos[1] > imgC.shape[1] / 2)
        pos[1] -= imgC.shape[1];
    phi r = {0, static_cast<double>(-pos[0]), static_cast<double>(-pos[1])};
    delete[] pos;
    return r;
}

phi rotateFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    arrayReal<double> corr = polarCorrelation(imgX, imgC);
    double angle = bestAngle(corr);
    phi r = {angle, 0, 0};
    return r;
}

/// version in python
imageComplex polarFFT(const imageReal<double> &img) {
    imageComplex f = fftShift2(fft2(img));
    imageReal<double> fr = f.abs();
    //    int *pos = fr.argmin();
    //    double r = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    //    delete[] pos;
    //    imageReal<int> mask = circularMask(img.shape, r);
    //    fr = fr * mask;
    for (int i = 0; i < 5; ++i)
        fr = hanning(fr);
    for (int i = 0; i < 1; ++i)
        fr = fr - hanning(fr);
    imageReal<double> p = polarBack(fr, true);
    return fft2(p);
}

phi alignFFT(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageReal<double> imgXD = denoiseSS(imgX);
    imageReal<double> imgCD = denoiseSS(imgC);
    imageComplex fx = polarFFT(imgXD);
    imageComplex fc = polarFFT(imgCD);
    imageComplex phaseF = fx * fc.conj() / (fx * fc.conj()).abs();
    imageReal<double> phase = ifft2C(phaseF).abs();
    int *pos = phase.argmax();
    double angle = static_cast<double>(pos[0]) / imgC.shape[0] * 360;
    delete[] pos;
    return angleAll(angle, imgX, imgC);
}

phi angleAll(double angle, const imageReal<double> &imgX, const imageReal<double> &imgC, bool fineFlag) {
    imageReal<double> imgC1 = rotate(imgC, angle);
    phi s1 = shiftFFTX(imgX, imgC1);
    imageReal<double> imgX1 = shift(rotate(imgC, angle), s1.x, s1.y);
    phi r = {angle, s1.x, s1.y};
    double corr1 = correntropy(imgX1, imgX);
    double corr = corr1;

    if (fineFlag) {
        for (int a = static_cast<int>(angle) - 20; a < static_cast<int>(angle) + 20; a += 2) {
            imageReal<double> imgC2 = rotate(imgC, a);
            phi s2 = shiftFFTX(imgX, imgC2);
            imageReal<double> imgX2 = shift(rotate(imgC, a), s2.x, s2.y);
            double corr2 = correntropy(imgX2, imgX);
            if (corr < corr2) {
                r = {static_cast<double>(a), s2.x, s2.y};
                corr = corr2;
            }
        }

        int angleP = static_cast<int>(r.angle);

        for (int a = angleP - 1; a < angleP + 2; ++a) {
            imageReal<double> imgC2 = rotate(imgC, a);
            phi s2 = shiftFFTX(imgX, imgC2);
            imageReal<double> imgX2 = shift(imgC2, s2.x, s2.y);
            double corr2 = correntropy(imgX2, imgX);
            if (corr < corr2) {
                r = {static_cast<double>(a), s2.x, s2.y};
                corr = corr2;
            }
        }
    }

    double angleT = r.angle > 180 ? r.angle - 180 : r.angle + 180;
    imageReal<double> imgC2 = rotate(imgC, angleT);
    phi s2 = shiftFFTX(imgX, imgC2);
    imageReal<double> imgX2 = shift(rotate(imgC, angleT), s2.x, s2.y);
    double corr2 = correntropy(imgX2, imgX);
    if (corr < corr2) {
        r = {angleT, s2.x, s2.y};
        corr = corr2;
    }

    if (std::abs(r.x) > 20 or std::abs(r.y) > 20)
        r.x = r.y = 0;
    return r;
}

phi alignFFT2(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageReal<double> fx = fftShift2(fft2(imgX)).abs();
    imageReal<double> fc = fftShift2(fft2(imgC)).abs();
    arrayReal<double> corr = polarCorrelation(fx, fc);
    double angle = bestAngle(corr);
    return angleAll(angle, imgX, imgC);
}

phi alignShiftPolar(const imageReal<double> &imgX, const imageReal<double> &imgC, int maxShiftX, int maxShiftY) {
    double corr = 0;
    phi r{};
    for (int i = 0; i < maxShiftX; ++i)
        for (int j = 0; j < maxShiftY; ++j) {
            imageReal<double> imgXN = shift(imgX, -i, -j);
            phi angle = rotateFFT(imgXN, imgC);
            imageReal<double> imgCN = shift(rotate(imgC, angle.angle), i, j);
            double corrT = correntropy(imgX, imgCN);
            if (corrT > corr) {
                corr = corrT;
                r = {angle.angle, static_cast<double>(i), static_cast<double>(j)};
            }
        }
    return r;
}

phi alignFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageReal<double> tempsr = imgC;
    phi final{};
    phi psr;
    for (int i = 0; i < 3; ++i) {
        phi pTempS, pTempR;
        pTempS = shiftFFTX(imgX, tempsr);
        tempsr = shift(tempsr, pTempS.x, pTempS.y);
        psr = psr + pTempS;
        pTempR = rotateFFTX(imgX, tempsr);
        tempsr = rotate(tempsr, pTempR.angle);
        psr = psr + pTempR;
        while (psr.angle > 360)
            psr.angle -= 360;
        if (std::abs(pTempS.x) < ALIGN_TH and std::abs(pTempS.y) < ALIGN_TH and pTempR.angle < ALIGN_TH)
            break;
    }
    imageReal<double> temprs = imgC;
    phi prs;
    for (int i = 0; i < 3; ++i) {
        phi pTempS, pTempR;
        pTempR = rotateFFTX(imgX, temprs);
        temprs = rotate(temprs, pTempR.angle);
        prs = prs + pTempR;
        pTempS = shiftFFTX(imgX, temprs);
        temprs = shift(temprs, pTempS.x, pTempS.y);
        prs = prs + pTempS;
        while (prs.angle > 360)
            prs.angle -= 360;
        if (std::abs(pTempS.x) < ALIGN_TH and std::abs(pTempS.y) < ALIGN_TH and pTempR.angle < ALIGN_TH)
            break;
    }
    imageReal<double> imgSR = applyPhi(imgC, psr);
    imageReal<double> imgRS = applyPhi(imgC, prs);
    if (correntropy(imgX, imgSR) > correntropy(imgX, imgRS))
        final = psr;
    else
        final = prs;
    if (std::abs(final.x) > 20 or std::abs(final.y) > 20)
        final.x = final.y = 0;
    return final;
}

imageReal<double> applyPhi(const imageReal<double> &img, phi p) {
    imageReal<double> r = rotate(img, p.angle);
    r = shift(r, p.x, p.y);
    return r;
}

imageReal<double> inversePhi(const imageReal<double> &img, phi p) {
    imageReal<double> r = shift(img, -p.x, -p.y);
    r = rotate(r, -p.angle);
    return r;
}

imageReal<double>
reference(const stackReal<double> &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &)) {
    if (stk.shape[0] == 1)
        return stk.pieceGet(0);
    int l = stk.shape[0] % 2 == 1 ? stk.shape[0] / 2 + 1 : stk.shape[0] / 2;
    int sp[3] = {l, stk.shape[1], stk.shape[2]};
    stackReal<double> s = stackReal<double>(sp);
    int count = 0;
    for (int i = 0; i < stk.shape[0]; i += 2) {
        imageReal<double> a = stk.pieceGet(i);
        if (i + 1 >= stk.shape[0]) {
            break;
        }
        imageReal<double> b = stk.pieceGet(i + 1);
        phi t = func(a, b);
        imageReal<double> bn = shift(rotate(b, t.angle), t.x, t.y);
        imageReal<double> avr = (a + bn) / 2;
        if (i % 10 == 0)
            avr = centerImage(avr);
        s.pieceSet(count, avr);
        ++count;
    }
    return reference(s, func);
}

imageReal<double>
reference(const stackReal<double> &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool),
          int width, bool fineFlag) {
    if (stk.shape[0] == 1)
        return stk.pieceGet(0);
    int l = stk.shape[0] % 2 == 1 ? stk.shape[0] / 2 + 1 : stk.shape[0] / 2;
    int sp[3] = {l, stk.shape[1], stk.shape[2]};
    stackReal<double> s = stackReal<double>(sp);
    int count = 0;
    for (int i = 0; i < stk.shape[0]; i += 2) {
        imageReal<double> a = stk.pieceGet(i);
        if (i + 1 >= stk.shape[0]) {
            s.pieceSet(count, a);
            ++count;
            break;
        }
        imageReal<double> b = stk.pieceGet(i + 1);
        phi t = func(a, b, width, fineFlag);
        imageReal<double> bn = shift(rotate(b, t.angle), t.x, t.y);
        imageReal<double> avr = (a + bn) / 2;
        if (i % 10 == 0)
            avr = centerImage(avr);
        s.pieceSet(count, avr);
        ++count;
    }
    return reference(s, func, width, fineFlag);
}

imageReal<double> reference(const stackReal<double> &stk,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double, bool),
                            int outWidth, int inWidth, double th, bool fineFlag) {
    if (stk.shape[0] == 1)
        return stk.pieceGet(0);
    int l = stk.shape[0] % 2 == 1 ? stk.shape[0] / 2 + 1 : stk.shape[0] / 2;
    int sp[3] = {l, stk.shape[1], stk.shape[2]};
    stackReal<double> s = stackReal<double>(sp);
    int count = 0;
    for (int i = 0; i < stk.shape[0]; i += 2) {
        imageReal<double> a = stk.pieceGet(i);
        if (i + 1 >= stk.shape[0]) {
            s.pieceSet(count, a);
            ++count;
            break;
        }
        imageReal<double> b = stk.pieceGet(i + 1);
        phi t = func(a, b, outWidth, inWidth, th, fineFlag);
        imageReal<double> bn = shift(rotate(b, t.angle), t.x, t.y);
        imageReal<double> avr = (a + bn) / 2;
        if (i % 10 == 0)
            avr = centerImage(avr);
        s.pieceSet(count, avr);
        ++count;
    }
    return reference(s, func, outWidth, inWidth, th, fineFlag);
}

imageReal<double> reference(const fileStack &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &)) {
    if (stk.shape[0] == 1)
        return stk.pieceGet(0);
    fileStack s = fileStack(stk.auxFolder, stk.auxFolderBack);
    int count = 0;
    for (int i = 0; i < stk.shape[0]; i += 2) {
        imageReal<double> a = stk.pieceGet(i);
        if (i + 1 >= stk.shape[0]) {
            break;
        }
        imageReal<double> b = stk.pieceGet(i + 1);
        phi t = func(a, b);
        imageReal<double> bn = shift(rotate(b, t.angle), t.x, t.y);
        imageReal<double> avr = (a + bn) / 2;
        if (i % 10 == 0)
            avr = centerImage(avr);
        s.append(avr);
        ++count;
    }
    if (access(s.auxFolder.c_str(), 0) == 0)
        system(("rm -r " + s.auxFolder).c_str());
    rename(s.auxFolderBack.c_str(), s.auxFolder.c_str());
    int mdB = mkdir(s.auxFolderBack.c_str(), S_IRWXU);
    if (mdB == -1)
        throw baseException("Error: Unable to create aux folder for files stacks!");
    return reference(s, func);
}

imageReal<double> centerImage(const imageReal<double> &img) {
    double m = img.mean();
    imageReal<double> imgN = img - m;
    imageReal<int> mask = circularMask(img.shape, img.shape[0] / 2.0);
    imageReal<double> imgC = imgN.copy();
    imageReal<double> imgX, imgY, imgXY;
    phi pFinal{};

    for (int iter = 0; iter < 5; ++iter) {
        imgC = imgC * mask;
        imgX = imgC.flipX();
        imgY = imgC.flipY();
        imgXY = imgC.flipAll();
        double shiftX = 0, shiftY = 0, nx = 0, ny = 0;
        phi p{};

        p = shiftFFTX(imgY, imgC);
        if (std::abs(p.y) < img.shape[1] / 3.0) {
            shiftY += p.y;
            ny++;
        }

        p = shiftFFTX(imgX, imgC);
        if (std::abs(p.x) < img.shape[0] / 3.0) {
            shiftX += p.x;
            nx++;
        }

        p = shiftFFTX(imgXY, imgC);
        if (std::abs(p.x) < img.shape[0] / 3.0) {
            shiftX += p.x;
            nx++;
        }
        if (std::abs(p.y) < img.shape[1] / 3.0) {
            shiftY += p.y;
            ny++;
        }

        shiftX = nx > 0 ? shiftX / nx : 0;
        shiftY = ny > 0 ? shiftY / ny : 0;
        pFinal.x += shiftX / 2;
        pFinal.y += shiftY / 2;
        imgC = applyPhi(imgN, pFinal);
        imgC = imgC * mask;

        p = rotateFFTX(imgY, imgC);
        while (p.angle > 360)
            p.angle -= 360;
        pFinal.angle = p.angle / 2;
        imgC = applyPhi(imgN, pFinal);
        imgC = imgC * mask;

        arrayReal<double> lineX = arrayReal<double>(imgC.shape[0]);
        arrayReal<double> lineY = arrayReal<double>(imgC.shape[1]);
        for (int i = 0; i < lineX.length; ++i)
            lineX[i] = imgC.column(i).max();
        for (int i = 0; i < lineY.length; ++i)
            lineY[i] = imgC.row(i).max();
        double thX = lineX.min() + 0.75 * (lineX.max() - lineX.min());
        double thY = lineY.min() + 0.75 * (lineY.max() - lineY.min());
        int x0 = 0, y0 = 0, xf = static_cast<int>(lineX.length) - 1, yf = static_cast<int>(lineY.length) - 1;
        while (lineX[x0] < thX)
            x0++;
        while (lineY[y0] < thY)
            y0++;
        while (lineX[xf] < thX)
            xf--;
        while (lineY[yf] < thY)
            yf--;
        if (xf - x0 > yf - y0) {
            pFinal.angle -= 90;
            imgC = applyPhi(imgN, pFinal);
        }
        if (std::abs(shiftX) < ALIGN_TH and std::abs(shiftY) < ALIGN_TH and p.angle < ALIGN_TH)
            break;
    }

    imgN = imgC;
    imgN = imgN + m;
    return imgN;
}

imageReal<double> initReference(const stackReal<double> &stk) {
    int shape[2] = {stk.shape[1], stk.shape[2]};
    imageReal<double> ref = imageReal<double>(shape);
    for (int i = 0; i < stk.shape[0]; ++i)
        ref = ref + stk.pieceGet(i);
    ref = ref / stk.shape[0];
    ref = centerImage(ref);

    return ref;
}

imageReal<double> initReference(const fileStack &stk) {
    int shape[2] = {stk.shape[1], stk.shape[2]};
    imageReal<double> ref = imageReal<double>(shape);
    for (int i = 0; i < stk.shape[0]; ++i)
        ref = ref + stk.pieceGet(i);
    ref = ref / stk.shape[0];
    ref = centerImage(ref);

    return ref;
}

void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &), bool updateFlag) {
    int N = stk.shape[0];
    double lambda = 1.0 / (N / 2.0);
    double lambdaP = 1 - lambda;
    phi p;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        p = func(ref, piece);
        imageReal<double> imgAl = applyPhi(piece, p);
        if (updateFlag) {
            ref = ref * lambdaP + imgAl * lambda;
            if (i % 10 == 0)
                ref = centerImage(ref);
        }
        stk.pieceSet(i, imgAl);
    }
}

void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool), int width, bool fineFlag,
                bool updateFlag) {
    int N = stk.shape[0];
    double lambda = 1.0 / (N / 2.0);
    double lambdaP = 1 - lambda;
    phi p;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        //        p = func(piece, ref, width, fineFlag);
        //        imageReal<double> imgAl = inversePhi(piece, p);
        p = func(ref, piece, width, fineFlag);
        imageReal<double> imgAl = applyPhi(piece, p);
        if (updateFlag) {
            ref = ref * lambdaP + imgAl * lambda;
            if (i % 10 == 0)
                ref = centerImage(ref);
        }
        stk.pieceSet(i, imgAl);
    }
}

void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double, bool), int outWidth,
                int inWidth, double th, bool fineFlag, bool updateFlag) {
    int N = stk.shape[0];
    double lambda = 1.0 / (N / 2.0);
    double lambdaP = 1 - lambda;
    phi p;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        //        p = func(piece, ref, outWidth, inWidth, th, fineFlag);
        //        imageReal<double> imgAl = inversePhi(piece, p);
        p = func(ref, piece, outWidth, inWidth, th, fineFlag);
        imageReal<double> imgAl = applyPhi(piece, p);
        if (updateFlag) {
            ref = ref * lambdaP + imgAl * lambda;
            if (i % 10 == 0)
                ref = centerImage(ref);
        }
        stk.pieceSet(i, imgAl);
    }
}

void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &), std::vector<std::string> info,
             const std::string &xmdName) {
    phi p;
    std::vector<std::string> xmd;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        p = func(ref, piece);
        imageReal<double> imgAl = applyPhi(piece, p);
        double corr = correntropy(imgAl, ref);
        std::string temp = info[i] + "     " + std::to_string(p.y) + "     " + std::to_string(p.x) + "     " +
                           std::to_string(p.angle) + "     " +
                           std::to_string(0) + "     " + std::to_string(corr) + ' ';
        xmd.emplace_back(temp);
    }
    writeXMD(xmd, xmdName);
}

void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool),
             std::vector<std::string> info, const std::string &xmdName, int width, bool fineFlag) {
    phi p;
    std::vector<std::string> xmd;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        p = func(ref, piece, width, fineFlag);
        imageReal<double> imgAl = applyPhi(piece, p);
        double corr = correntropy(imgAl, ref);
        std::string temp = info[i] + "     " + std::to_string(p.y) + "     " + std::to_string(p.x) + "     " +
                           std::to_string(p.angle) + "     " +
                           std::to_string(0) + "     " + std::to_string(corr) + ' ';
        xmd.emplace_back(temp);
    }
    writeXMD(xmd, xmdName);
}

void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double th, bool fineFlag),
             std::vector<std::string> info, const std::string &xmdName, int outWidth, int inWidth, double th,
             bool fineFlag) {
    phi p;
    std::vector<std::string> xmd;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        p = func(ref, piece, outWidth, inWidth, th, fineFlag);
        imageReal<double> imgAl = applyPhi(piece, p);
        double corr = correntropy(imgAl, ref);
        std::string temp = info[i] + "     " + std::to_string(p.y) + "     " + std::to_string(p.x) + "     " +
                           std::to_string(p.angle) + "     " +
                           std::to_string(0) + "     " + std::to_string(corr) + ' ';
        xmd.emplace_back(temp);
    }
    writeXMD(xmd, xmdName);
}

void
refinement(fileStack &stk, imageReal<double> &ref, phi (*func)(const imageReal<double> &, const imageReal<double> &)) {
    int N = stk.shape[0];
    double lambda = 1.0 / (N / 2.0);
    double lambdaP = 1 - lambda;
    phi p;
    for (int i = 0; i < stk.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        p = func(ref, piece);
        imageReal<double> imgAl = applyPhi(piece, p);
        ref = ref * lambdaP + imgAl * lambda;
        if (i % 10 == 0)
            ref = centerImage(ref);
    }
}

void refinementWithAux(stackReal<double> &stk, stackReal<double> &stkAux, imageReal<double> &ref) {
    int N = stk.shape[0];
    double lambda = 1.0 / (N / 2.0);
    double lambdaP = 1 - lambda;
    phi p;
    for (int i = 0; i < stkAux.shape[0]; ++i) {
        imageReal<double> piece = stk.pieceGet(i);
        imageReal<double> pieceAux = stkAux.pieceGet(i);
        p = alignFFTX(ref, pieceAux);
        imageReal<double> imgAl = applyPhi(piece, p);
        imageReal<double> imgAlAux = applyPhi(pieceAux, p);
        ref = ref * lambdaP + imgAlAux * lambda;
        if (i % 10 == 0)
            ref = centerImage(ref);
        stk.pieceSet(i, imgAl);
        stkAux.pieceSet(i, imgAlAux);
    }
}

phi alignPower(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageComplex fx = fftShift2(fft2(imgX));
    imageComplex fc = fftShift2(fft2(imgC));
    imageReal<double> spx = fx.abs();
    spx = spx.flipX();
    imageReal<double> spc = fc.abs();
    int centerX = imgX.shape[0] / 2;
    int centerY = imgX.shape[1] / 2;
    imageReal<double> delta = spx * spx - spc * spc;
    imageReal<int> deltaZero = imageReal<int>(imgX.shape);
    arrayReal<int> angleGroup = arrayReal<int>(360);
    double value = delta.abs().max();
    for (int i = 0; i < delta.shape[0]; ++i)
        for (int j = 0; j < delta.shape[1]; ++j) {
            if (std::abs(delta.get(i, j)) > value / 10) {
                deltaZero.set(i, j, 1);
                int angle = static_cast<int>(std::round(std::atan2(centerX - i, j - centerY) / M_PI / 2 * 360));
                angle = angle < 0 ? angle + 180 : angle;
                angleGroup[angle]++;
            }
        }
    int pos = static_cast<int>(angleGroup.argmax());
    phi p{};
    p.angle = pos;
    return p;
}

phi gravityCenter(const imageReal<double> &img, int x, int y, double rate = 0.8) {
    bool flag = true;
    int width = 1;
    double s = 0, sX = 0, sY = 0;
    double valueMax = img.get(x, y);
    while (flag) {
        for (int i = x - width; i <= x + width; ++i)
            for (int j = y - width; j <= y + width; ++j) {
                double value = img.get(i, j);
                if (0 <= i and i < img.shape[0] and 0 <= j and j < img.shape[1]) {
                    if (value > valueMax * rate) {
                        s += value;
                        sX += value * i;
                        sY += value * j;
                    } else {
                        flag = false;
                    }
                }
            }
        ++width;
    }
    double xF = sX / s;
    double yF = sY / s;
    phi r = {
            0,
            xF,
            yF,
    };
    return r;
}

phi alignPeak(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageReal<double> fx = fftShift2(fft2(imgX)).abs();
    imageReal<double> fc = fftShift2(fft2(imgC)).abs();
    int centerX = fx.shape[0] / 2;
    int centerY = fx.shape[1] / 2;
    fx = fx * fx / (fx.get(centerX, centerY) * fx.get(centerX, centerY));
    fc = fc * fc / (fc.get(centerX, centerY) * fc.get(centerX, centerY));
    int width = imgX.shape[0] / 20;
    for (int i = centerX - width; i <= centerX + width; ++i)
        for (int j = centerY - width; j <= centerY + width; ++j) {
            fx.set(i, j, 0);
            fc.set(i, j, 0);
        }
    int *posX = fx.argmax();
    int xX = posX[0];
    int yX = posX[1];
    delete[] posX;
    int *posC = fc.argmax();
    int xC = posC[0];
    int yC = posC[1];
    delete[] posC;
    double angleX = std::atan2(centerX - xX, yX - centerY) / 2 / M_PI * 360;
    angleX = angleX < 0 ? angleX + 360 : angleX;
    double angleC = std::atan2(centerX - xC, yC - centerY) / 2 / M_PI * 360;
    angleC = angleC < 0 ? angleC + 360 : angleC;
    double angleD = angleX - angleC < 0 ? angleX - angleC + 360 : angleX - angleC;
    return angleAll(angleD, imgX, imgC);
}

int circleEstimatePeak(const stackReal<double> &stk, bool fineFlag) {
    imageReal<double> ref = stk.pieceGet(0);
    int r = 0;
    double v = 0;
    int start = stk.shape[1] / 40;
    int end = stk.shape[1] / 10;
    for (int i = start; i < end; ++i) {
        imageReal<double> refTemp = ref.copy();
        for (int j = 1; j < stk.shape[0]; ++j) {
            imageReal<double> img = stk.pieceGet(j);
            phi p = alignPeak(ref, img, i, fineFlag);
            imageReal<double> imgAl = applyPhi(img, p);
            refTemp = refTemp + imgAl;
        }
        refTemp = refTemp / stk.shape[0];
        double vT = correntropy(refTemp, ref);
        if (vT > v) {
            v = vT;
            r = i;
        }
    }
    return r;
}

int circleEstimateShape(const stackReal<double> &stk, int inWidth, double th, bool fineFlag) {
    imageReal<double> ref = stk.pieceGet(0);
    int r = 0;
    double v = 0;
    int start = stk.shape[1] / 40;
    int end = stk.shape[1] / 10;
    for (int i = start; i < end; ++i) {
        imageReal<double> refTemp = ref.copy();
        for (int j = 1; j < stk.shape[0]; ++j) {
            imageReal<double> img = stk.pieceGet(j);
            phi p = alignShape(ref, img, i, inWidth, th, fineFlag);
            imageReal<double> imgAl = applyPhi(img, p);
            refTemp = refTemp + imgAl;
        }
        refTemp = refTemp / stk.shape[0];
        double vT = correntropy(refTemp, ref);
        if (vT > v) {
            v = vT;
            r = i;
        }
    }
    return r;
}

/** Gravity center is not used as default. You may uncomment code in this function and comment the
 * the corresponding part to use gravity center.
 *
 * Note the variable width is the parameter corresponding to high pass filter in power spectrum domain
 * which could be estimated using function circleEstimate or be manually estimated.
*/
phi alignPeak(const imageReal<double> &imgX, const imageReal<double> &imgC, int width, bool fineFlag) {
    imageReal<double> fx = fftShift2(fft2(imgX)).abs();
    imageReal<double> fc = fftShift2(fft2(imgC)).abs();
    int centerX = fx.shape[0] / 2;
    int centerY = fx.shape[1] / 2;
    fx = fx * fx / (fx.get(centerX, centerY) * fx.get(centerX, centerY));
    fc = fc * fc / (fc.get(centerX, centerY) * fc.get(centerX, centerY));
    if (width == -1)
        width = imgX.shape[0] / 15;
    for (int i = centerX - width; i <= centerX + width; ++i)
        for (int j = centerY - width; j <= centerY + width; ++j) {
            if ((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY) <= width * width) {
                fx.set(i, j, 0);
                fc.set(i, j, 0);
            }
        }
    int *posX = fx.argmax();
    int xX = posX[0];
    int yX = posX[1];
    //    phi pX = gravityCenter(fx, xX, yX);
    //    double xXG = pX.x;
    //    double yXG = pX.y;
    delete[] posX;
    int *posC = fc.argmax();
    int xC = posC[0];
    int yC = posC[1];
    //    phi pC = gravityCenter(fc, xC, yC);
    //    double xCG = pC.x;
    //    double yCG = pC.y;
    delete[] posC;
    //    double angleX = std::atan2(centerX - xXG, yXG - centerY) / 2 / M_PI * 360;
    double angleX = std::atan2(centerX - xX, yX - centerY) / 2 / M_PI * 360;
    angleX = angleX < 0 ? angleX + 360 : angleX;
    //    double angleC = std::atan2(centerX - xCG, yCG - centerY) / 2 / M_PI * 360;
    double angleC = std::atan2(centerX - xC, yC - centerY) / 2 / M_PI * 360;
    angleC = angleC < 0 ? angleC + 360 : angleC;
    double angleD = angleX - angleC < 0 ? angleX - angleC + 360 : angleX - angleC;
    return angleAll(angleD, imgX, imgC, fineFlag);
}

phi alignShape(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    imageReal<double> fx = fftShift2(fft2(imgX)).abs();
    imageReal<double> fc = fftShift2(fft2(imgC)).abs();
    int centerX = fx.shape[0] / 2;
    int centerY = fx.shape[1] / 2;
    fx = fx * fx * fx;
    fc = fc * fc * fc;
    //    int outWidth = fx.shape[0] / 15;
    int outWidth = 22;
    int inWidth = 0;
    int cornerSize = fx.shape[0] / 10;
    imageReal<double> cornerX = fx.get(0, 0, cornerSize, cornerSize);
    imageReal<double> cornerC = fc.get(0, 0, cornerSize, cornerSize);
    double rate = 1;
    double vx = cornerX.mean() * rate;
    double vc = cornerC.mean() * rate;
    for (int i = 0; i < fx.shape[0]; ++i)
        for (int j = 0; j < fx.shape[1]; ++j) {
            int distance = (i - centerX) * (i - centerX) + (j - centerY) * (j - centerY);
            if (distance <= inWidth * inWidth or distance >= outWidth * outWidth or fx.get(i, j) < vx)
                fx.set(i, j, 0);
            if (distance <= inWidth * inWidth or distance >= outWidth * outWidth or fc.get(i, j) < vc)
                fc.set(i, j, 0);
        }
    double angleX = momentDirection(fx, outWidth);
    double angleC = momentDirection(fc, outWidth);
    double angleD = angleX - angleC < 0 ? angleX - angleC + 360 : angleX - angleC;
    return angleAll(angleD, imgX, imgC);
}

phi alignShape(const imageReal<double> &imgX, const imageReal<double> &imgC, int outWidth, int inWidth, double th,
               bool fineFlag) {
    imageReal<double> fx = fftShift2(fft2(imgX)).abs();
    imageReal<double> fc = fftShift2(fft2(imgC)).abs();
    int centerX = fx.shape[0] / 2;
    int centerY = fx.shape[1] / 2;
    fx = fx * fx * fx;
    fc = fc * fc * fc;
    if (inWidth == -1)
        inWidth = 0;
    double vx, vc;
    if (th == -1) {
        int cornerSize = fx.shape[0] / 10;
        double rate = 1;
        imageReal<double> cornerX = fx.get(0, 0, cornerSize, cornerSize);
        imageReal<double> cornerC = fc.get(0, 0, cornerSize, cornerSize);
        vx = cornerX.mean() * rate;
        vc = cornerC.mean() * rate;
    } else
        vx = vc = th;

    for (int i = 0; i < fx.shape[0]; ++i)
        for (int j = 0; j < fx.shape[1]; ++j) {
            int distance = (i - centerX) * (i - centerX) + (j - centerY) * (j - centerY);
            if (distance <= inWidth * inWidth or distance >= outWidth * outWidth or fx.get(i, j) < vx)
                fx.set(i, j, 0);
            if (distance <= inWidth * inWidth or distance >= outWidth * outWidth or fc.get(i, j) < vc)
                fc.set(i, j, 0);
        }
    double angleX = momentDirection(fx, outWidth);
    double angleC = momentDirection(fc, outWidth);
    double angleD = angleX - angleC < 0 ? angleX - angleC + 360 : angleX - angleC;
    return angleAll(angleD, imgX, imgC, fineFlag);
}

template<typename T>
imageReal<int> thresholdAdjust(const imageReal<T> &img, T value) {
    imageReal<int> r = imageReal<int>(img.shape);
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j) {
            if (img.get(i, j) < value)
                r.set(i, j, 0);
            else
                r.set(i, j, 1);
        }
    return r;
}

phi alignParticle(const imageReal<double> &imgX, const imageReal<double> &imgC) {
    int size = 320;
    arrayReal<int> recX = coordParticle(imgX, size);
    arrayReal<int> recC = coordParticle(imgC, size);
    double xX = static_cast<double>(recX[0] + recX[2]) / 2;
    double yX = static_cast<double>(recX[1] + recX[3]) / 2;
    double xC = static_cast<double>(recC[0] + recC[2]) / 2;
    double yC = static_cast<double>(recC[1] + recC[3]) / 2;
    double d0X = std::sqrt((xX - recX[0]) * (xX - recX[0]) + (yX - recX[1]) * (yX - recX[1]));
    double d0C = std::sqrt((xC - recC[0]) * (xC - recC[0]) + (yC - recC[1]) * (yC - recC[1]));
    double centerX[2] = {xX, yX};
    double centerC[2] = {xC, yC};
    imageReal<double> px = butterworthLow(imgX, centerX, d0X, 10);
    imageReal<double> pc = butterworthLow(imgC, centerC, d0C, 10);
    imageReal<double> fx = fftShift2(fft2(px)).abs();
    imageReal<double> fc = fftShift2(fft2(pc)).abs();
    double mx = fx.mean();
    double mc = fc.mean();
    imageReal<int> fxb = thresholdAdjust(fx, 4 * mx);
    imageReal<int> fcb = thresholdAdjust(fc, 4 * mc);
    fxb = largestDomain(fxb);
    fcb = largestDomain(fcb);
    double angleX = momentDirection(fxb);
    angleX = angleX < 0 ? angleX + 360 : angleX;
    double angleC = momentDirection(fcb);
    angleC = angleC < 0 ? angleC + 360 : angleC;
    double angleD = angleX - angleC < 0 ? angleX - angleC + 360 : angleX - angleC;
    return angleAll(angleD, imgX, imgC);
}
