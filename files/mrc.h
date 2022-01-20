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

#ifndef EM_MRC_H
#define EM_MRC_H

/** @file
 * this file contains classes and functions for reading/writing MRC file
 * reference website:
 *      https://bio3d.colorado.edu/imod/doc/mrc_format.txt
 *      https://www.ccpem.ac.uk/mrc_format/mrc2014.php
*/

#include <fstream>
#include <cstring>
#include <iostream>
#include "stack3D.h"

/// flag for different data types in MRC files
#define MODE_INT_8 0
#define MODE_INT_16 1
#define MODE_REAL_32 2
#define MODE_COMPLEX_16 3
#define MODE_COMPLEX_32 4
#define MODE_UINT_16 6

/** @brief
 * class for processing head information in MRC files
 * see variable meaning in reference website
*/
class mrcHeader {
private:
    int totalSize = 1024; // header information takes total size of 1024 bytes

public:
    int nx;
    int ny;
    int nz;
    int mode;
    int nxStart;
    int nyStart;
    int nzStart;
    int mx;
    int my;
    int mz;
    int cellA[3]{};
    int cellB[3]{};
    int mapC;
    int mapR;
    int mapS;
    double dMin;
    double dMax;
    double dMean;
    int ispg;
    int nSymBT;
    int extra[25]{};
    char extTyp[4];
    int nVersion;
    int origin[3]{};
    char map[4];
    unsigned char machSt[4];
    int rms;
    int nLabl;
    int label[200]{};

public:
    mrcHeader();

    mrcHeader(const mrcHeader &header);

    void read(std::ifstream &mrc);

    void write(std::ofstream &mrc) const;
};

/** @brief
 * class for processing overall information in MRC files
 * doesn't consider datatype conversion from big endian machines to little endian machines or vice versa
 * currently only image stack is supported
 * TODO: add support to volume data
*/
class mrcFile {
public:
    mrcHeader header; // header information
    stackReal<double> data; // imageStack if MRC file is real data type
    stackComplex dataComplex; // imageStack if MRC file is complex data type
    bool real; // whether MRC file is real data type

private:
    /// change dimension of ary according to mid, like transpose in pytorch
    template<typename T>
    void arrayArrange(T *ary, const int *mid, int len) {
        T *temp = new T[len];
        for (int i = 0; i < len; ++i)
            temp[i] = ary[i];
        for (int i = 0; i < len; ++i) {
            ary[i] = temp[mid[i]];
        }
        delete[] temp;
    }

    /// read data from file stream
    void read(std::ifstream &mrc);

    /// write data to file stream
    void write(std::ofstream &mrc);

public:
    mrcFile();

    mrcFile(const mrcFile &mrc);

    explicit mrcFile(const stackReal<double> &stk);

    explicit mrcFile(const imageReal<double> &img);

    explicit mrcFile(const stackComplex &stk);

    explicit mrcFile(const imageComplex &img);

    /// read data from file path
    void read(const std::string &fileName);

    /// update header information according to data
    void update();

    /// read data to file path
    void write(const std::string &fileName);

    /// set real image stack stk to data
    void setData(const stackReal<double> &stk);

    /// set real complex stack stk to data
    void setData(const stackComplex &stk);
};

#endif //EM_MRC_H
