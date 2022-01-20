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

#include "xmippXMD.h"

stackReal<double> readXMD(const std::string &path, std::vector<std::string> &info) {
    std::ifstream in;
    in.open(path, std::ios::in);
    if (!in)
        throw baseException("Error: Unable to open txt files!");
    std::string line;
    std::vector<std::vector<std::string>> r;
    int row = 0;
    while (getline(in, line)) {
        if (row < 5) {
            ++row;
            continue;
        }
        std::istringstream sin(line);
        std::vector<std::string> fields;
        std::string field;
        while (getline(sin, field, '@')) {
            fields.emplace_back(field);
        }
        if (fields.size() != 2)
            break;
        r.emplace_back(fields);
        info.emplace_back(fields[0] + '@' + fields[1]);
    }

    in.close();
    int sliceNumber = static_cast<int>(r.size());
    std::set<std::string> stkNames;
    for (int i = 0; i < sliceNumber; ++i) {
        std::vector<std::string> piece = r[i];
        stkNames.insert(piece[1]);
    }

    stackReal<double> stks[stkNames.size()];

    std::set<std::string>::iterator iter;
    int ii;
    for (ii = 0, iter = stkNames.begin(); ii < stkNames.size(); ++ii, ++iter) {
        mrcFile mrc = mrcFile();
        std::string name = *iter;
        if (name[name.length() - 1] == ' ')
            name = name.substr(0, name.length() - 1);
        mrc.read(name);
        stks[ii] = mrc.data;
    }

    imageReal<double> imgStk[sliceNumber];
    for (int i = 0; i < sliceNumber; ++i) {
        std::vector<std::string> piece = r[i];
        int index = atoi(piece[0].c_str()) - 1;
        std::string stkName = piece[1];
        for (ii = 0, iter = stkNames.begin(); ii < stkNames.size(); ++ii, ++iter) {
            if (stkName == *iter) {
                imgStk[i] = stks[ii].pieceGet(index);
                break;
            }
        }
    }

    stackReal<double> final = image2Stack(imgStk, sliceNumber);
    return final;
}

std::string zfill(const std::string &number, int length) {
    std::string ss = number;
    while (ss.size() < length) {
        ss.insert(0, "0");
    }
    return ss;
}

void generateInfo(const std::string &path, int number, std::vector<std::string> &info) {
    for (int i = 1; i <= number; ++i) {
        std::string index = zfill(std::to_string(i), 6);
        info.emplace_back(index.append("@" + path));
    }
}

void writeXMD(const std::vector<std::string> &info, const std::string &name) {
    std::string path = name + "_alignment.xmd";
    std::ofstream fout(path, std::ios::out);
    fout << "# XMIPP_STAR_1 * " << std::endl;
    fout << "# " << std::endl;
    fout << "data_noname" << std::endl;
    fout << "loop_" << std::endl;
    fout << " _image" << std::endl;
    fout << " _shiftX" << std::endl;
    fout << " _shiftY" << std::endl;
    fout << " _anglePsi" << std::endl;
    fout << " _flip" << std::endl;
    fout << " _maxCC" << std::endl;
    for (const auto &line: info) {
        fout << line << std::endl;
    }
    fout.close();
}