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

#include "fft.h"

arrayComplex fft(const arrayReal<double> &ary)
{
    int N = static_cast<int>(ary.length);
    double *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<double *>(fftw_malloc(sizeof(double) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
        memcpy(in, ary.data, N * sizeof(double));

        fftw_execute(p);
        fftw_destroy_plan(p);

        int index = N / 2 + 1;

        for (int i = index; i < N; ++i)
        {
            out[i][0] = out[N - i][0];
            out[i][1] = -out[N - i][1];
        }
    }

    arrayComplex aryOut = arrayComplex(out, N);

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

arrayComplex fft(const arrayComplex &ary)
{
    int N = static_cast<int>(ary.length);
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        for (int i = 0; i < N; ++i)
        {
            in[i][0] = ary.data[2 * i];
            in[i][1] = ary.data[2 * i + 1];
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, N);

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

arrayComplex ifftC(const arrayComplex &ary)
{
    int N = static_cast<int>(ary.length);
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (int i = 0; i < N; ++i)
        {
            in[i][0] = ary.data[2 * i];
            in[i][1] = ary.data[2 * i + 1];
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, N) / N;

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

arrayReal<double> ifftR(const arrayComplex &ary)
{
    int N = static_cast<int>(ary.length);
    fftw_complex *in;
    int index = N / 2 + 1;
    double *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * index));
    out = static_cast<double *>(fftw_malloc(sizeof(double) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);

        for (int i = 0; i < index; ++i)
        {
            in[i][0] = ary.data[2 * i];
            in[i][1] = ary.data[2 * i + 1];
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayReal<double> aryOut = arrayReal<double>(out, N) / N;

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

arrayComplex fftShift(const arrayComplex &ary)
{
    arrayComplex aryOut = arrayComplex(ary.length);
    int half0 = static_cast<int>(ary.length) / 2;
    int half1 = ary.length % 2 == 1 ? half0 + 1 : half0;
    memcpy(aryOut.data, ary.data + 2 * half1, 2 * half0 * sizeof(double));
    memcpy(aryOut.data + 2 * half0, ary.data, 2 * half1 * sizeof(double));
    return aryOut;
}

arrayComplex ifftShift(const arrayComplex &ary)
{
    arrayComplex aryOut = arrayComplex(ary.length);
    int half0 = static_cast<int>(ary.length) / 2;
    int half1 = ary.length % 2 == 1 ? half0 + 1 : half0;
    memcpy(aryOut.data, ary.data + 2 * half0, 2 * half1 * sizeof(double));
    memcpy(aryOut.data + 2 * half1, ary.data, 2 * half0 * sizeof(double));
    return aryOut;
}

arrayComplex fftHalf(const arrayReal<double> &ary)
{
    int N = static_cast<int>(ary.length);
    int width = N / 2 + 1;
    double *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<double *>(fftw_malloc(sizeof(double) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * width));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
        memcpy(in, ary.data, N * sizeof(double));

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, width) / N;

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

arrayReal<double> ifftHalf(const arrayComplex &ary, int originLength)
{
    int N = static_cast<int>(ary.length);
    int width = originLength % 2 == 1 ? 2 * N - 1 : 2 * N - 2;
    fftw_complex *in;
    double *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<double *>(fftw_malloc(sizeof(double) * width));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_c2r_1d(width, in, out, FFTW_ESTIMATE);

        for (int i = 0; i < N; ++i)
        {
            in[i][0] = ary.data[2 * i];
            in[i][1] = ary.data[2 * i + 1];
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayReal<double> aryOut = arrayReal<double>(out, width);

    fftw_free(in);
    fftw_free(out);
    return aryOut;
}

imageComplex fft2(const imageReal<double> &img)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    int width = img.shape[1] / 2 + 1;
    size_t N2 = static_cast<size_t>(img.shape[0]) * width;
    double *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<double *>(fftw_malloc(sizeof(double) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N2));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_r2c_2d(img.shape[0], img.shape[1], in, out, FFTW_ESTIMATE);

        memcpy(in, img.data.data, N * sizeof(double));

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    fftw_complex *outN;
    outN = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    if (outN == nullptr)
        throw baseException("Error: Insufficient available memory for FFT!");

    for (int i = 0; i < img.shape[0]; ++i)
        for (int j = 0; j < img.shape[1]; ++j)
        {
            if (j < width)
            {
                outN[i * img.shape[1] + j][0] = out[i * (img.shape[1] / 2 + 1) + j][0];
                outN[i * img.shape[1] + j][1] = out[i * (img.shape[1] / 2 + 1) + j][1];
            }
            else if (i == 0)
            {
                outN[i * img.shape[1] + j][0] = out[i * (img.shape[1] / 2 + 1) + img.shape[1] - j][0];
                outN[i * img.shape[1] + j][1] = -out[i * (img.shape[1] / 2 + 1) + img.shape[1] - j][1];
            }
            else
            {
                outN[i * img.shape[1] + j][0] = out[(img.shape[0] - i) * (img.shape[1] / 2 + 1) + img.shape[1] - j][0];
                outN[i * img.shape[1] + j][1] = -out[(img.shape[0] - i) * (img.shape[1] / 2 + 1) + img.shape[1] - j][1];
            }
        }

    arrayComplex aryOut = arrayComplex(outN, N);
    imageComplex imgOut = imageComplex(aryOut, img.shape);

    fftw_free(in);
    fftw_free(out);
    fftw_free(outN);

    return imgOut;
}

imageComplex fft2(const imageComplex &img)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_2d(img.shape[0], img.shape[1], in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        for (size_t i = 0; i < N; ++i)
        {
            double *temp = img.data[i];
            in[i][0] = temp[0];
            in[i][1] = temp[1];
            delete[] temp;
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, N);
    imageComplex imgOut = imageComplex(aryOut, img.shape);

    fftw_free(in);
    fftw_free(out);
    return imgOut;
}

imageComplex ifft2C(const imageComplex &img)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_2d(img.shape[0], img.shape[1], in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (size_t i = 0; i < N; ++i)
        {
            double *temp = img.data[i];
            in[i][0] = temp[0];
            in[i][1] = temp[1];
            delete[] temp;
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, N);
    imageComplex imgOut = imageComplex(aryOut, img.shape) / static_cast<double>(N);

    fftw_free(in);
    fftw_free(out);
    return imgOut;
}

imageReal<double> ifft2R(const imageComplex &img)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    int width = img.shape[1] / 2 + 1;
    size_t N2 = static_cast<size_t>(img.shape[0]) * width;
    fftw_complex *in;
    double *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N2));
    out = static_cast<double *>(fftw_malloc(sizeof(double) * N));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_c2r_2d(img.shape[0], img.shape[1], in, out, FFTW_ESTIMATE);

        int count = 0;
        for (int i = 0; i < img.shape[0]; i++)
            for (int j = 0; j < width; j++)
            {
                double *value = img.get(i, j);
                in[count][0] = value[0];
                in[count][1] = value[1];
                delete[] value;
                count++;
            }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayReal<double> aryOut = arrayReal<double>(out, N);
    imageReal<double> imgOut = imageReal<double>(aryOut, img.shape) / N;

    fftw_free(in);
    fftw_free(out);
    return imgOut;
}

imageComplex fftShift2(const imageComplex &img)
{
    imageComplex imgOut = imageComplex(img.shape);
    int d0 = img.shape[0] % 2 == 1 ? img.shape[0] / 2 + 1 : img.shape[0] / 2;
    int d1 = img.shape[1] % 2 == 1 ? img.shape[1] / 2 + 1 : img.shape[1] / 2;
    imageComplex img00 = img.get(0, 0, d0, d1);
    imageComplex img01 = img.get(0, d1, d0, img.shape[1] - d1);
    imageComplex img10 = img.get(d0, 0, img.shape[0] - d0, d1);
    imageComplex img11 = img.get(d0, d1, img.shape[0] - d0, img.shape[1] - d1);
    imgOut.set(img.shape[0] - d0, img.shape[1] - d1, img00);
    imgOut.set(img.shape[0] - d0, 0, img01);
    imgOut.set(0, img.shape[1] - d1, img10);
    imgOut.set(0, 0, img11);
    return imgOut;
}

imageComplex ifftShift2(const imageComplex &img)
{
    imageComplex imgOut = imageComplex(img.shape);
    int d0 = img.shape[0] % 2 == 1 ? img.shape[0] / 2 + 1 : img.shape[0] / 2;
    int d1 = img.shape[1] % 2 == 1 ? img.shape[1] / 2 + 1 : img.shape[1] / 2;
    imageComplex img00 = img.get(0, 0, img.shape[0] - d0, img.shape[1] - d1);
    imageComplex img01 = img.get(img.shape[0] - d0, 0, d0, img.shape[1] - d1);
    imageComplex img10 = img.get(0, img.shape[1] - d1, img.shape[0] - d0, d1);
    imageComplex img11 = img.get(img.shape[0] - d0, img.shape[1] - d1, d0, d1);
    imgOut.set(d0, d1, img00);
    imgOut.set(0, d1, img01);
    imgOut.set(d0, 0, img10);
    imgOut.set(0, 0, img11);
    return imgOut;
}

imageComplex fft2Half(const imageReal<double> &img)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    int width = img.shape[1] / 2 + 1;
    int shapeOut[2] = {img.shape[0], width};
    size_t N2 = static_cast<size_t>(shapeOut[0]) * shapeOut[1];
    double *in;
    fftw_complex *out;
    fftw_plan p;

    in = static_cast<double *>(fftw_malloc(sizeof(double) * N));
    out = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N2));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_r2c_2d(img.shape[0], img.shape[1], in, out, FFTW_ESTIMATE);

        memcpy(in, img.data.data, N * sizeof(double));

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayComplex aryOut = arrayComplex(out, N2) / static_cast<double>(N);
    imageComplex imgOut = imageComplex(aryOut, shapeOut);

    fftw_free(in);
    fftw_free(out);
    return imgOut;
}

imageReal<double> ifft2Half(const imageComplex &img, int originLength)
{
    size_t N = static_cast<size_t>(img.shape[0]) * img.shape[1];
    int width = originLength % 2 == 1 ? 2 * img.shape[1] - 1 : 2 * (img.shape[1] - 1);
    int shapeOut[2] = {img.shape[0], width};
    size_t N2 = static_cast<size_t>(shapeOut[0]) * shapeOut[1];
    fftw_complex *in;
    double *out;
    fftw_plan p;

    in = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N));
    out = static_cast<double *>(fftw_malloc(sizeof(double) * N2));

    if (in == nullptr or out == nullptr)
    {
        throw baseException("Error: Insufficient available memory for FFT!");
    }
    else
    {
        p = fftw_plan_dft_c2r_2d(shapeOut[0], shapeOut[1], in, out, FFTW_ESTIMATE);
        for (size_t i = 0; i < N; ++i)
        {
            double *temp = img.data[i];
            in[i][0] = temp[0];
            in[i][1] = temp[1];
            delete[] temp;
        }

        fftw_execute(p);
        fftw_destroy_plan(p);
    }

    arrayReal<double> aryOut = arrayReal<double>(out, N2);
    imageReal<double> imgOut = imageReal<double>(aryOut, shapeOut);

    fftw_free(in);
    fftw_free(out);
    return imgOut;
}