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

#ifndef ALIGNMENT_XMIPPXMD_H
#define ALIGNMENT_XMIPPXMD_H


#include <fstream>
#include <sstream>
#include <set>
#include "stack3D.h"
#include "mrcFile.h"


stackReal<float> readXMD(const std::string &path, std::vector<std::string> &info);


std::string zfill(const std::string &number, int length);


void generateInfo(const std::string &path, int number, std::vector<std::string> &info);


void writeXMD(const std::vector<std::string> &info, const std::string &name);


#endif //ALIGNMENT_XMIPPXMD_H