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

#ifndef ALIGNMENT_ALIGNMENTEVALUATION_H
#define ALIGNMENT_ALIGNMENTEVALUATION_H

#include "alignFFT.h"
#include "mrcFile.h"
#include "matrix.h"
#include <fstream>
#include <cstring>
#include <sstream>

/// Data in MRC format
stackReal<float> readData(const std::string &path);

/// Alignment parameters stored in txt file, in format "index angle shiftX shiftY"
std::vector<std::vector<float>> readParams(const std::string &path);

void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout);

void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout);

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath);

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath);

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout);

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout);

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath);

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath);

void evaluateRefImage(const std::string &dataPath, const std::string &refPath, std::ofstream &fout);

#endif // ALIGNMENT_ALIGNMENTEVALUATION_H