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

#ifndef ALIGNMENT_NOISE_H
#define ALIGNMENT_NOISE_H

#include "image2D.h"
#include "fft.h"

float randomReal(float minValue, float maxValue);

int randomInt(int minValue, int maxValue);

template <typename T>
imageReal<float> gaussianNoise(const imageReal<T> &img, double rate)
{
    double m1 = img.min();
    double m2 = img.max();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    imageReal<float> r;
    if (typeid(T) != typeid(float))
        r = img.copy().template asType<float>();
    else
        r = img.copy();
    float var = rate * img.var();

    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0, std::sqrt(var));
    for (int i = 0; i < r.data.length; ++i)
    {
        auto gr = static_cast<float>(distribution(generator));
        float value = r.data[i] + gr;
        r.data[i] = value > m2 ? m2 : value;
        r.data[i] = value < m1 ? m1 : r.data[i];
    }
    return r;
}

template <typename T>
imageReal<float> denoiseSS(const imageReal<T> &img)
{
    imageComplex f = fftShift2(fft2(img));
    imageReal<double> fa = f.abs();
    imageReal<double> ang = f.angle();
    fa = fa * fa;
    int width = fa.shape[0] / 10;
    double value = 0;
    value += fa.get(0, 0, width, width).mean();
    value += fa.get(fa.shape[0] - width, 0, width, width).mean();
    value += fa.get(0, fa.shape[1] - width, width, width).mean();
    value += fa.get(fa.shape[0] - width, fa.shape[1] - width, width, width).mean();
    value /= 4;
    fa = fa - value;
    fa = fa.clip(0, INFINITY);
    fa = fa.sqrt();
    imageComplex nf = bind(fa, ang);
    nf = ifftShift2(nf);
    imageReal<float> r = ifft2C(nf).abs().asType<float>();
    return r;
}

template <typename T1, typename T2>
imageReal<float> convolutionSpace(const imageReal<T1> &img, const imageReal<T2> &kernel)
{
    if (kernel.shape[0] % 2 == 0 || kernel.shape[1] % 2 == 0)
        throw baseException("Error: Size of convolution kernel should be odd integer!");
    imageReal<float> r = imageReal<float>(img.shape);
    int h = (kernel.shape[0] - 1) / 2, w = (kernel.shape[1] - 1) / 2;
    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            float value = 0;
            for (int m = -h; m <= h; ++m)
                for (int n = -w; n <= w; ++n)
                {
                    int kx = h + m, ky = w + n, ix = i + m, iy = j + n;
                    if (ix < 0)
                        ix += img.shape[0];
                    else if (ix > img.shape[0] - 1)
                        ix -= img.shape[0];
                    if (iy < 0)
                        iy += img.shape[1];
                    else if (iy > img.shape[1] - 1)
                        iy -= img.shape[1];
                    value += img.get(ix, iy) * kernel.get(kx, ky);
                }
            r.set(i, j, value);
        }
    return r;
}

template <typename T1, typename T2>
imageReal<float> convolutionFourier(const imageReal<T1> &img, const imageReal<T2> &kernel)
{
    if (kernel.shape[0] % 2 == 0 || kernel.shape[1] % 2 == 0)
        throw baseException("Error: Size of convolution kernel should be odd integer!");
    int h = (kernel.shape[0] - 1) / 2, w = (kernel.shape[1] - 1) / 2;
    imageReal<T1> imgP = padding(img, h, h, w, w, PAD_WARP);
    //    imgP = padding(imgP, 0, imgP.shape[0] - 2 * h, 0, imgP.shape[1] - 2 * w, PAD_CONSTANT);
    imgP = padding(imgP, 0, h, 0, w, PAD_CONSTANT);
    imageReal<T2> kernelP = imageReal<T2>(imgP.shape);
    for (int i = h; i < 2 * h + 1; ++i)
        for (int j = w; j < 2 * w + 1; ++j)
            kernelP.set(i - h, j - w, kernel.get(i, j));
    for (int i = h; i < 2 * h + 1; ++i)
        for (int j = 0; j < w; ++j)
            kernelP.set(i - h, kernelP.shape[1] - w + j, kernel.get(i, j));
    for (int i = 0; i < h; ++i)
        for (int j = w; j < 2 * w + 1; ++j)
            kernelP.set(kernelP.shape[0] - h + i, j - w, kernel.get(i, j));
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < w; ++j)
            kernelP.set(kernelP.shape[0] - h + i, kernelP.shape[1] - w + j, kernel.get(i, j));

    imageComplex fi = fft2(imgP);
    imageComplex fk = fft2(kernelP);
    imageComplex fr = fi * fk.conj();
    imageReal<double> rP = ifft2C(fr).real();

    rP = rP.get(h, w, img.shape[0], img.shape[1]);

    imageReal<float> r = rP.asType<float>();

    return r;
}

template <typename T>
imageReal<float> gaussianBlur(const imageReal<T> &img, int size = 3)
{
    double sigma = 0.3 * ((size - 1) * 0.5 - 1) + 0.8;
    int s[2] = {size, size};
    imageReal<double> kernel = imageReal<double>(s);
    int c = (size + 1) / 2;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            double value = 1 / (2 * M_PI * sigma * sigma) * exp(-((i - c) * (i - c) + (j - c) * (j - c)) / (2 * sigma * sigma));
            kernel.set(i, j, value);
        }
    kernel = kernel / kernel.sum();
    return convolutionFourier(img, kernel);
}

#endif //ALIGNMENT_NOISE_H