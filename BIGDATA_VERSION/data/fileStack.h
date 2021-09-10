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

#ifndef ALIGNMENT_FILESTACK_H
#define ALIGNMENT_FILESTACK_H


#include "mrcFile.h"
#include "stack3D.h"
#include <vector>
#include <glob.h>
#include <sstream>
#include "image2D.h"
#include <zconf.h>
#include <sys/stat.h>


std::vector<std::string> split(const std::string &line, char symbol);


class fileStack {
public:
    std::vector<std::string> names;
    std::string auxFolder;
    std::string auxFolderBack;
    int shape[3]{};

    fileStack() = default;

    explicit fileStack(const std::string &aux, const std::string &auxB);

    fileStack(const std::string &path, const std::string &extension, const std::string &aux,
              const std::string &auxB);

    float get(int i, int j, int k) const;

    imageReal<float> pieceGet(int index) const;

    void pieceSet(int index, const imageReal<float> &img);

    void append(const std::string &path);

    void append(const imageReal<float> &img);


};


#endif //ALIGNMENT_FILESTACK_H