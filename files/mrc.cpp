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

#include "mrc.h"

mrcHeader::mrcHeader() {
    nx = 0;
    ny = 0;
    nz = 0;
    mode = MODE_REAL_32;
    nxStart = 0;
    nyStart = 0;
    nzStart = 0;
    mx = 0;
    my = 0;
    mz = 0;
    for (int &i: cellA) {
        i = 0;
    }
    for (int &i: cellB) {
        i = 0;
    }
    mapC = 1;
    mapR = 2;
    mapS = 3;
    dMin = 0;
    dMax = 0;
    dMean = 0;
    ispg = 0;
    nSymBT = 0;
    for (int &i: extra) {
        i = 0;
    }
    for (char &i: extTyp) {
        i = '\000';
    }
    nVersion = 0;
    for (int &i: origin) {
        i = 0;
    }
    map[0] = 'M';
    map[1] = 'A';
    map[2] = 'P';
    map[3] = ' ';
    machSt[0] = 0x44;
    machSt[1] = 0x44;
    machSt[2] = 0x00;
    machSt[3] = 0x00;
    rms = 0;
    nLabl = 0;
    for (int &i: label) {
        i = 0;
    }
}

mrcHeader::mrcHeader(const mrcHeader &header) {
    nx = header.nx;
    ny = header.ny;
    nz = header.nz;
    mode = header.mode;
    nxStart = header.nxStart;
    nyStart = header.nyStart;
    nzStart = header.nzStart;
    mx = header.mx;
    my = header.my;
    mz = header.mz;
    for (int i = 0; i < 3; ++i)
        cellA[i] = header.cellA[i];
    for (int i = 0; i < 3; ++i)
        cellB[i] = header.cellB[i];
    mapC = header.mapC;
    mapR = header.mapR;
    mapS = header.mapS;
    dMin = header.dMin;
    dMax = header.dMax;
    dMean = header.dMean;
    ispg = header.ispg;
    nSymBT = header.nSymBT;
    for (int i = 0; i < 25; ++i)
        extra[i] = header.extra[i];
    for (int i = 0; i < 4; ++i)
        extTyp[i] = header.extTyp[i];
    nVersion = header.nVersion;
    for (int i = 0; i < 3; ++i)
        origin[i] = header.origin[i];
    for (int i = 0; i < 4; ++i)
        map[i] = header.map[i];
    for (int i = 0; i < 4; ++i)
        machSt[i] = header.machSt[i];
    rms = header.rms;
    nLabl = header.nLabl;
    for (int i = 0; i < 200; ++i)
        label[i] = header.label[i];
}

void mrcHeader::read(std::ifstream &mrc) {
    char headerData[totalSize];
    mrc.read(headerData, totalSize);
    int *head = reinterpret_cast<int *>(headerData);
    nx = *(head + 0);
    ny = *(head + 1);
    nz = *(head + 2);
    mode = *(head + 3);
    nxStart = *(head + 4);
    nyStart = *(head + 5);
    nzStart = *(head + 6);
    mx = *(head + 7);
    my = *(head + 8);
    mz = *(head + 9);
    memcpy(cellA, head + 10, 3 * sizeof(int));
    memcpy(cellB, head + 13, 3 * sizeof(int));
    mapC = *(head + 16);
    mapR = *(head + 17);
    mapS = *(head + 18);
    dMin = static_cast<double>(*(reinterpret_cast<float *>(headerData) + 19));
    dMax = static_cast<double>(*(reinterpret_cast<float *>(headerData) + 20));
    dMean = static_cast<double>(*(reinterpret_cast<float *>(headerData) + 21));
    ispg = *(head + 22);
    nSymBT = *(head + 23);
    memcpy(extra, head + 24, 25 * sizeof(int));
    memcpy(extTyp, headerData + 104, 4);
    nVersion = *(head + 27);
    memcpy(origin, head + 49, 3 * sizeof(int));
    memcpy(map, headerData + 208, 4);
    memcpy(machSt, headerData + 212, 4);
    rms = *(head + 54);
    nLabl = *(head + 55);
    memcpy(label, head + 56, 200 * sizeof(int));
}

void mrcHeader::write(std::ofstream &mrc) const {
    mrc.write((char *) &nx, sizeof(int));
    mrc.write((char *) &ny, sizeof(int));
    mrc.write((char *) &nz, sizeof(int));
    int modeOut = MODE_REAL_32;
    mrc.write((char *) &modeOut, sizeof(int));
    mrc.write((char *) &nxStart, sizeof(int));
    mrc.write((char *) &nyStart, sizeof(int));
    mrc.write((char *) &nzStart, sizeof(int));
    mrc.write((char *) &mx, sizeof(int));
    mrc.write((char *) &my, sizeof(int));
    mrc.write((char *) &mz, sizeof(int));
    mrc.write((char *) &cellA, 3 * sizeof(int));
    mrc.write((char *) &cellB, 3 * sizeof(int));
    mrc.write((char *) &mapC, sizeof(int));
    mrc.write((char *) &mapR, sizeof(int));
    mrc.write((char *) &mapS, sizeof(int));
    mrc.write((char *) &dMin, sizeof(float));
    mrc.write((char *) &dMax, sizeof(float));
    mrc.write((char *) &dMean, sizeof(float));
    mrc.write((char *) &ispg, sizeof(int));
    mrc.write((char *) &nSymBT, sizeof(int));
    mrc.write((char *) &extra, 25 * sizeof(int));
    mrc.write((char *) &origin, 3 * sizeof(int));
    mrc.write(map, 4);
    mrc.write((char *) &machSt, sizeof(int));
    mrc.write((char *) &rms, sizeof(int));
    mrc.write((char *) &nLabl, sizeof(int));
    mrc.write((char *) &label, 200 * sizeof(int));
}

void mrcFile::read(std::ifstream &mrc) {
    header.read(mrc);
    int axis[3];
    axis[0] = header.mapC;
    axis[1] = header.mapR;
    axis[2] = header.mapS;
    if (header.map[0] != 'M' or header.map[1] != 'A' or header.map[2] != 'P' or header.map[3] != ' ')
        throw baseException("Error: MAP identifier of files type is wrong!");
    if (header.mapC != 1 and header.mapR != 2 and header.mapS != 3) {
        std::cout << "Note: File not saved in Column-Row-Stack order and it will be converted to CRS order!"
                  << std::endl;
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
    if (header.mode == MODE_INT_8 or header.mode == MODE_INT_16 or header.mode == MODE_REAL_32 or
        header.mode == MODE_UINT_16) {
        real = true;
        data = stackReal<double>(dataShape);
    } else if (header.mode == MODE_COMPLEX_16 or header.mode == MODE_COMPLEX_32) {
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
                        data.set(i, j, k, static_cast<double>(*value));
                        break;
                    }
                    case MODE_INT_16: {
                        auto *value = reinterpret_cast<int16_t *>(temp);
                        data.set(i, j, k, static_cast<double>(*value));
                        break;
                    }
                    case MODE_REAL_32: {
                        auto *value = reinterpret_cast<float *>(temp);
                        data.set(i, j, k, *value);
                        break;
                    }
                    case MODE_UINT_16: {
                        auto *value = reinterpret_cast<uint16_t *>(temp);
                        data.set(i, j, k, static_cast<double>(*value));
                        break;
                    }
                    case MODE_COMPLEX_16: {
                        auto *valueR = reinterpret_cast<int16_t *>(temp);
                        mrc.read(temp, step);
                        auto *valueI = reinterpret_cast<int16_t *>(temp);
                        double value[2] = {static_cast<double>(*valueR), static_cast<double>(*valueI)};
                        dataComplex.set(i, j, k, value);
                    }
                    case MODE_COMPLEX_32: {
                        auto *valueR = reinterpret_cast<float *>(temp);
                        mrc.read(temp, step);
                        auto *valueI = reinterpret_cast<float *>(temp);
                        double value[2] = {static_cast<double>(*valueR), static_cast<double>(*valueI)};
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
                    auto value = static_cast<float>(data.get(i, j, k));
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
    data = stackReal<double>();
    dataComplex = stackComplex();
    real = true;
}

mrcFile::mrcFile(const mrcFile &mrc) {
    header = mrcHeader(mrc.header);
    data = stackReal<double>(mrc.data);
    dataComplex = stackComplex(mrc.dataComplex);
    real = mrc.real;
}

mrcFile::mrcFile(const stackReal<double> &stk) {
    header = mrcHeader();
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<double>(stk);
    dataComplex = stackComplex();
    real = true;
    update();
}

mrcFile::mrcFile(const imageReal<double> &img) {
    header = mrcHeader();
    stackReal<double> stk = image2Stack(img);
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<double>(stk);
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
    data = stackReal<double>();
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
    data = stackReal<double>();
    dataComplex = stackComplex(stk);
    real = false;
    update();
}

void mrcFile::read(const std::string &fileName) {
    std::ifstream mrc;
    mrc.open(fileName, std::ios::binary | std::ios::in);
    if (!mrc)
        throw baseException("Error: Unable to read MRC files of " + fileName + "!");
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
        throw baseException("Error: Unable to write MRC files of " + fileName + "!");
    write(mrc);
    mrc.close();
}

void mrcFile::setData(const stackReal<double> &stk) {
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<double>(stk);
    real = true;
    dataComplex = stackComplex();
    update();
}

void mrcFile::setData(const stackComplex &stk) {
    header.nx = stk.shape[2];
    header.ny = stk.shape[1];
    header.nz = stk.shape[0];
    data = stackReal<double>();
    real = false;
    dataComplex = stackComplex(stk);
    update();
}
