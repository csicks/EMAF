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

#include "mrcHeader.h"

mrcHeader::mrcHeader()
{
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
    for (int &i : cellA)
    {
        i = 0;
    }
    for (int &i : cellB)
    {
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
    for (int &i : extra)
    {
        i = 0;
    }
    for (char &i : extTyp)
    {
        i = '\000';
    }
    nVersion = 0;
    for (int &i : origin)
    {
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
    for (int &i : label)
    {
        i = 0;
    }
}

mrcHeader::mrcHeader(const mrcHeader &header)
{
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

void mrcHeader::read(std::ifstream &mrc)
{
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
    dMin = *(reinterpret_cast<float *>(headerData) + 19);
    dMax = *(reinterpret_cast<float *>(headerData) + 20);
    dMean = *(reinterpret_cast<float *>(headerData) + 21);
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

void mrcHeader::write(std::ofstream &mrc) const
{
    mrc.write((char *)&nx, sizeof(int));
    mrc.write((char *)&ny, sizeof(int));
    mrc.write((char *)&nz, sizeof(int));
    int modeOut = MODE_REAL_32;
    mrc.write((char *)&modeOut, sizeof(int));
    mrc.write((char *)&nxStart, sizeof(int));
    mrc.write((char *)&nyStart, sizeof(int));
    mrc.write((char *)&nzStart, sizeof(int));
    mrc.write((char *)&mx, sizeof(int));
    mrc.write((char *)&my, sizeof(int));
    mrc.write((char *)&mz, sizeof(int));
    mrc.write((char *)&cellA, 3 * sizeof(int));
    mrc.write((char *)&cellB, 3 * sizeof(int));
    mrc.write((char *)&mapC, sizeof(int));
    mrc.write((char *)&mapR, sizeof(int));
    mrc.write((char *)&mapS, sizeof(int));
    mrc.write((char *)&dMin, sizeof(float));
    mrc.write((char *)&dMax, sizeof(float));
    mrc.write((char *)&dMean, sizeof(float));
    mrc.write((char *)&ispg, sizeof(int));
    mrc.write((char *)&nSymBT, sizeof(int));
    mrc.write((char *)&extra, 25 * sizeof(int));
    mrc.write((char *)&origin, 3 * sizeof(int));
    mrc.write(map, 4);
    mrc.write((char *)&machSt, sizeof(int));
    mrc.write((char *)&rms, sizeof(int));
    mrc.write((char *)&nLabl, sizeof(int));
    mrc.write((char *)&label, 200 * sizeof(int));
}