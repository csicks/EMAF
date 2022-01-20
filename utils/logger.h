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

#ifndef EM_LOGGER_H
#define EM_LOGGER_H

/** @file
 * this file contains classes and functions for logger and exceptions
*/

#include <fstream>
#include <iostream>

/** @brief
 * class for exceptions
 * @example throw baseException("Error");
*/
class baseException : std::exception {
public:
    std::string msg; // string to print

public:
    explicit baseException(std::string s);
};

/** @brief
 * class for logger which outputs to both screen and file
*/
class logger {
public:
    std::string name; // path to record log

public:
    explicit logger(std::string s);

    /// write string to both screen and file
    void write(const std::string &s) const;
};

#endif //EM_LOGGER_H
