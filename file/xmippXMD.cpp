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

#include "xmippXMD.h"

stackReal<float> readXMD(const std::string &path, std::vector<std::string> &info)
{
    std::ifstream in;
    in.open(path, std::ios::in);
    if (!in)
        throw baseException("Error: Unable to open txt file!");
    std::string line;
    std::vector<std::vector<std::string>> r;
    int row = 0;
    while (getline(in, line))
    {
        if (row < 5)
        {
            ++row;
            continue;
        }
        std::istringstream sin(line);
        std::vector<std::string> fields;
        std::string field;
        while (getline(sin, field, '@'))
        {
            fields.emplace_back(field);
        }
        if (fields.size() != 2)
            break;
        r.emplace_back(fields);
        info.emplace_back(fields[0] + '@' + fields[1]);
    }

    in.close();
    int sliceNumber = r.size();
    std::set<std::string> stkNames;
    for (int i = 0; i < sliceNumber; ++i)
    {
        std::vector<std::string> piece = r[i];
        stkNames.insert(piece[1]);
    }

    stackReal<float> stks[stkNames.size()];

    std::set<std::string>::iterator iter;
    int ii;
    for (ii = 0, iter = stkNames.begin(); ii < stkNames.size(); ++ii, ++iter)
    {
        mrcFile mrc = mrcFile();
        std::string name = *iter;
        if (name[name.length() - 1] == ' ')
            name = name.substr(0, name.length() - 1);
        mrc.read(name);
        stks[ii] = mrc.data;
    }

    imageReal<float> imgStk[sliceNumber];
    for (int i = 0; i < sliceNumber; ++i)
    {
        std::vector<std::string> piece = r[i];
        int index = atoi(piece[0].c_str()) - 1;
        std::string stkName = piece[1];
        for (ii = 0, iter = stkNames.begin(); ii < stkNames.size(); ++ii, ++iter)
        {
            if (stkName == *iter)
            {
                imgStk[i] = stks[ii].pieceGet(index);
                break;
            }
        }
    }

    stackReal<float> final = image2Stack(imgStk, sliceNumber);
    return final;
}

std::string zfill(const std::string &number, int length)
{
    std::string ss = number;
    while (ss.size() < length)
    {
        ss = "0" + ss;
    }
    return ss;
}

void generateInfo(const std::string &path, int number, std::vector<std::string> &info)
{
    for (int i = 1; i <= number; ++i)
    {
        std::string index = zfill(std::to_string(i), 6);
        info.emplace_back(index + '@' + path);
    }
}

void writeXMD(const std::vector<std::string> &info, const std::string &name)
{
    std::string path = name + "_alignment.xmd";
    std::ofstream outFile(path, std::ios::out);
    outFile << "# XMIPP_STAR_1 * " << std::endl;
    outFile << "# " << std::endl;
    outFile << "data_noname" << std::endl;
    outFile << "loop_" << std::endl;
    outFile << " _image" << std::endl;
    outFile << " _shiftX" << std::endl;
    outFile << " _shiftY" << std::endl;
    outFile << " _anglePsi" << std::endl;
    outFile << " _flip" << std::endl;
    outFile << " _maxCC" << std::endl;
    for (int i = 0; i < info.size(); ++i)
    {
        std::string line = info[i];
        outFile << line << std::endl;
    }
    outFile.close();
}