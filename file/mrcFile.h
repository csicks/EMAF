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

#ifndef ALIGNMENT_MRCFILE_H
#define ALIGNMENT_MRCFILE_H


#include <iostream>
#include "mrcHeader.h"
#include "stack3D.h"


/// Does not consider datatype conversion from big endian machines to little endian machines or vice versa
class mrcFile {
public:
    mrcHeader header;
    stackReal<float> data;
    stackComplex dataComplex;
    bool real;

private:
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

    void read(std::ifstream &mrc);

    void write(std::ofstream &mrc);

public:
    mrcFile();

    mrcFile(const mrcFile &mrc);

    explicit mrcFile(const stackReal<float> &stk);

    explicit mrcFile(const imageReal<float> &img);

    explicit mrcFile(const stackComplex &stk);

    explicit mrcFile(const imageComplex &img);

    void read(const std::string &fileName);

    void update();

    void write(const std::string &fileName);

    void setData(const stackReal<float> &stk);

    void setData(const stackComplex &stk);

};


#endif //ALIGNMENT_MRCFILE_H