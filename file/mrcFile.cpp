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

#include "mrcFile.h"


void mrcFile::read(std::ifstream &mrc) {
    header.read(mrc);
    int axis[3];
    axis[0] = header.mapC;
    axis[1] = header.mapR;
    axis[2] = header.mapS;
    if (header.map[0] != 'M' || header.map[1] != 'A' || header.map[2] != 'P' || header.map[3] != ' ')
        throw baseException("Error: MAP identifier of file type is wrong!");
    if (header.mapC != 1 && header.mapR != 2 && header.mapS != 3) {
        std::cout << "Note: File not saved in Column-Row-Stack order and it will be converted to CRS order!" << std::endl;
        int nxyz[3] = {header.nx, header.ny, header.nz};
        arrayArrange(nxyz, axis, 3);
        header.nx = nxyz[0];
        header.ny = nxyz[1];
        header.nz = nxyz[2];
        int nxyzStart[3] = {header.nxStart, header.nyStart, header.nzStart};
        arrayArrange(nxyzStart, axis, 3);
        header.nxStart = nxyzStart[0];
        header.nyStart = nxyzStart[1];
        header.nzStart = nxyzStart[2];
        int mxyz[3] = {header.mx, header.my, header.mz};
        arrayArrange(mxyz, axis, 3);
        header.mx = mxyz[0];
        header.my = mxyz[1];
        header.mz = mxyz[2];
        arrayArrange(header.cellA, axis, 3);
        arrayArrange(header.cellB, axis, 3);
    }

    int dataShape[3] = {header.nz, header.ny, header.nx};
    if (header.mode == MODE_INT_8 || header.mode == MODE_INT_16 || header.mode == MODE_REAL_32 || header.mode == MODE_UINT_16) {
        real = true;
        data = stackReal<float>(dataShape);
    } else if (header.mode == MODE_COMPLEX_16 || header.mode == MODE_COMPLEX_32) {
        real = false;
        dataComplex = stackComplex(dataShape);
    } else
        throw baseException("Error: MRC Mode value is invalid!");

    int step;
    switch (header.mode) {
        case MODE_INT_8:
            step = sizeof(char);
            break;
        case MODE_INT_16:
            step = 2 * sizeof(char);
            break;
        case MODE_REAL_32:
            step = 4 * sizeof(char);
            break;
        case MODE_COMPLEX_16:
            step = 2 * sizeof(char);
            break;
        case MODE_COMPLEX_32:
            step = 4 * sizeof(char);
            break;
        case MODE_UINT_16:
            step = 2 * sizeof(char);
            break;
    }

    char *temp = new char[step];
    for (int i = 0; i < data.shape[0]; ++i)
        for (int j = 0; j < data.shape[1]; ++j)
            for (int k = 0; k < data.shape[2]; ++k) {
                mrc.read(temp, step);
                switch (header.mode) {
                    case MODE_INT_8: {
                        auto *value = reinterpret_cast<int8_t *>(temp);
                        data.set(i, j, k, static_cast<float>(*value));
                        break;
                    }
                    case MODE_INT_16: {
                        auto *value = reinterpret_cast<int16_t *>(temp);
                        data.set(i, j, k, static_cast<float>(*value));
                        break;
                    }
                    case MODE_REAL_32: {
                        auto *value = reinterpret_cast<float *>(temp);
                        data.set(i, j, k, *value);
                        break;
                    }
                    case MODE_UINT_16: {
                        auto *value = reinterpret_cast<uint16_t *>(temp);
                        data.set(i, j, k, static_cast<float>(*value));
                        break;
                    }
                    case MODE_COMPLEX_16: {
                        auto *valueR = reinterpret_cast<int16_t *>(temp);
                        mrc.read(temp, step);
                        auto *valueI = reinterpret_cast<int16_t *>(temp);
                        double value[2] = {static_cast<double >(*valueR), static_cast<double >(*valueI)};
                        dataComplex.set(i, j, k, value);
                    }
                    case MODE_COMPLEX_32: {
                        auto *valueR = reinterpret_cast<float *>(temp);
                        mrc.read(temp, step);
                        auto *valueI = reinterpret_cast<float *>(temp);
                        double value[2] = {static_cast<double >(*valueR), static_cast<double >(*valueI)};
                        dataComplex.set(i, j, k, value);
                    }
                    default:
                        throw baseException("Error: Complex value is not supported now!");
                }
            }
    delete[] temp;
}


void mrcFile::write(std::ofstream &mrc) {
    update();
    header.write(mrc);
    for (int i = 0; i < data.shape[0]; ++i)
        for (int j = 0; j < data.shape[1]; ++j)
            for (int k = 0; k < data.shape[2]; ++k) {
                if (real) {
                    float value = data.get(i, j, k);
                    mrc.write((char *) &(value), sizeof(float));
                } else {
                    double *value = dataComplex.get(i, j, k);
                    auto valueR = static_cast<float>(value[0]);
                    auto valueI = static_cast<float>(value[1]);
                    delete[] value;
                    mrc.write((char *) &(valueR), sizeof(float));
                    mrc.write((char *) &(valueI), sizeof(float));
                }
            }
}


mrcFile::mrcFile() {
    header = mrcHeader();
    data = stackReal<float>();
    dataComplex = stackComplex();
    real = true;
}


mrcFile::mrcFile(const mrcFile &mrc) {
    header = mrcHeader(mrc.header);
    data = stackReal<float>(mrc.data);
    dataComplex = stackComplex(mrc.dataComplex);
    real = mrc.real;
}


mrcFile::mrcFile(const stackReal<float> &stk) {
    header = mrcHeader();
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>(stk);
    dataComplex = stackComplex();
    real = true;
    update();
}


mrcFile::mrcFile(const imageReal<float> &img) {
    header = mrcHeader();
    stackReal<float> stk = image2Stack(img);
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>(stk);
    dataComplex = stackComplex();
    real = true;
    update();
}


mrcFile::mrcFile(const stackComplex &stk) {
    header = mrcHeader();
    header.mode = MODE_COMPLEX_32;
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>();
    dataComplex = stackComplex(stk);
    real = false;
    update();
}


mrcFile::mrcFile(const imageComplex &img) {
    header = mrcHeader();
    header.mode = MODE_COMPLEX_32;
    stackComplex stk = image2Stack(img);
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>();
    dataComplex = stackComplex(stk);
    real = false;
    update();
}


void mrcFile::read(const std::string &fileName) {
    std::ifstream mrc;
    mrc.open(fileName, std::ios::binary | std::ios::in);
    if (!mrc)
        throw baseException("Error: Unable to read MRC file of " + fileName + "!");
    read(mrc);
    mrc.close();
    update();
}


void mrcFile::update() {
    if (real) {
        header.dMin = data.min();
        header.dMax = data.max();
        header.dMean = data.mean();
    } else {
        header.dMin = -1;
        header.dMax = -2;
        header.dMean = -3;
        header.rms = -1;
    }
}


void mrcFile::write(const std::string &fileName) {
    std::ofstream mrc;
    mrc.open(fileName, std::ios::binary | std::ios::out);
    if (!mrc)
        throw baseException("Error: Unable to write MRC file of " + fileName + "!");
    write(mrc);
    mrc.close();
}


void mrcFile::setData(const stackReal<float> &stk) {
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>(stk);
    real = true;
    dataComplex = stackComplex();
    update();
}


void mrcFile::setData(const stackComplex &stk) {
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<float>();
    real = false;
    dataComplex = stackComplex(stk);
    update();
}