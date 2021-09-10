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

/// Alignment parameters stored in txt file, in format 'index angle shiftX shiftY'
std::vector<std::vector<float>> readParams(const std::string &path);

/** For functions below, they are used for synthetic datasets results evaluation, with 'dataPath' correspondind to particle 
 *  stacks which is named as 'stack.mrcs' in each folder in synthetic datasets, 'txtPath' correspondind to transformation 
 *  recording text file which is named as 'parameters.txt' in each folder in synthetic datasets, and 'refPath' corresponding
 *  to origin image which is named 'ori.mrcs' in each synthetic dataset.
 * **/

/// Average reference error
void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout);

/// Average reference error
void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout);

/// Reference error for each image (for plotting error histogram)
void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath);

/// Reference error for each image (for plotting error histogram)
void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath);

/// Average pair-wise error (computed from nearby images in original order)
void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout);

/// Average pair-wise error (computed from nearby images in original order)
void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout);

/// Pair-wise error for each image (for plotting error histogram)
void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath);

/// Pair-wise error for each image (for plotting error histogram)
void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath);

/** The below function is used for real datasets evaluation, with 'dataPath' corresponding to real world particle stacks and
 *  'refPath' corresponding to the generated reference image.
 * **/
/// Average reference similarity
void evaluateRefImage(const std::string &dataPath, const std::string &refPath, std::ofstream &fout);

#endif //ALIGNMENT_ALIGNMENTEVALUATION_H