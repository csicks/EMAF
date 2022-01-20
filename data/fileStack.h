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

#ifndef EM_FILESTACK_H
#define EM_FILESTACK_H

/** @file
 * this file contains data structure for three-dimensional image stack which is stored in files instead of memory
 * only real value supported
*/

#include "mrc.h"
#include "stack3D.h"
#include <vector>
#include <glob.h>
#include <sstream>
#include "image2D.h"
#include <zconf.h>
#include <sys/stat.h>

/// split string
std::vector<std::string> split(const std::string &line, char symbol);

/** @brief
 * class for image stack which is stored in files, whose size is changeable
*/
class fileStack {
public:
    std::vector<std::string> names; // file names of images
    std::string auxFolder; // path of auxiliary folder 1 for file stack
    std::string auxFolderBack; // path of auxiliary folder 2 for file stack
    int shape[3]{}; // shape of file stack

    fileStack() = default;

    explicit fileStack(const std::string &aux, const std::string &auxB);

    /** construct from file stacks
     * @param path: path to input folder
     * @param extension: input file extension
    */
    fileStack(const std::string &path, const std::string &extension, const std::string &aux,
              const std::string &auxB);

    /// @return value at position (i,j,k)
    double get(int i, int j, int k) const;

    /// @return the ith piece of file stack
    imageReal<double> pieceGet(int index) const;

    /// set the ith piece of file stack to img
    void pieceSet(int index, const imageReal<double> &img);

    /// append image of path to the end of file stack
    void append(const std::string &path);

    /// append image to the end of file stack
    void append(const imageReal<double> &img);
};

#endif //EM_FILESTACK_H