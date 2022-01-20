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

#include "logger.h"

baseException::baseException(std::string s) : msg(std::move(s)) {
    std::cout << msg << std::endl;
}

logger::logger(std::string s) : name(std::move(s)) {
    std::ifstream f(name.c_str());
    if (f.good()) {
        int flag = remove(name.c_str());
        if (flag != 0)
            std::cout << "Error: fail to initialize the logger!" << std::endl;
    }
}

void logger::write(const std::string &s) const {
    std::ofstream fout(name, std::ios::out | std::ios::app);
    fout << s << std::endl;
    fout.close();
    std::cout << s << std::endl;
}