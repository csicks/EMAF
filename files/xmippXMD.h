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

#ifndef EM_XMIPPXMD_H
#define EM_XMIPPXMD_H

/** @file
 * this file contains functions for reading/writing XMD file in XMIPP
 * currently, only limited typed of XMD files are supported (those used in SPREAD)
*/

#include <fstream>
#include <sstream>
#include <set>
#include "stack3D.h"
#include "mrc.h"

/// read XMD file to get data and header information
stackReal<double> readXMD(const std::string &path, std::vector<std::string> &info);

/** fill string numbers with fixed length, like zfill in python
 * @example "12".zfill(5) -> "00012"
*/
std::string zfill(const std::string &number, int length);

/// generate information for certain XMD files
void generateInfo(const std::string &path, int number, std::vector<std::string> &info);

/// write data and header information to XMD file
void writeXMD(const std::vector<std::string> &info, const std::string &name);

#endif //EM_XMIPPXMD_H