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

#ifndef EM_ALIGNEVALUATION_H
#define EM_ALIGNEVALUATION_H

/** @file
 * this file mainly contains functions relevant to 2D image alignment method evaluation
*/

#include <fstream>
#include <cstring>
#include <sstream>
#include "align.h"
#include "mrc.h"
#include "matrix.h"

/** common used variable explanation:
 * @param dataPath: path for input image stacks
 * @param txtPath: path for ground truth alignment parameters
 * @param refPath: path for reference image if needed
 * @param func: given align method [function pointer]
 * @param width: width for high-pass/low-pass filter
 * @param outPath: path to store outputs
 * @param fout: output stream for files, which is an alternative choice when outPath specified before
*/

/// read image stack from given path, for data in MRC format
stackReal<double> readData(const std::string &path);

/// read alignment parameters stored in txt files, in format "index angle shiftX shiftY"
std::vector<std::vector<double>> readParams(const std::string &path);

/** evaluate error when aligned to the reference image
 * align function format: align(imgX, imgC) [suitable for all align methods]
 * @return average error
*/
void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi (*func)(const imageReal<double> &, const imageReal<double> &), std::ofstream &fout);

/** evaluate error when aligned to the reference image
 * align function format: align(imgX, imgC, width) [suitable for alignPeak and alignShape]
 * @return average error
*/
void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi (*func)(const imageReal<double> &, const imageReal<double> &, int), int width,
                       std::ofstream &fout);

/** evaluate error when aligned to the reference image
 * align function format: align(imgX, imgC) [suitable for all align methods]
 * @return error for every image
*/
void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi (*func)(const imageReal<double> &, const imageReal<double> &),
                           const std::string &outPath);

/** evaluate error when aligned to the reference image
 * align function format: align(imgX, imgC, width) [suitable for alignPeak and alignShape]
 * @return error for every image
*/
void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                           int width, const std::string &outPath);

/** evaluate error when aligned to neighborhood image (every two neighbor images to be aligned)
 * align function format: align(imgX, imgC) [suitable for all align methods]
 * @return average error
*/
void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi (*func)(const imageReal<double> &, const imageReal<double> &), std::ofstream &fout);

/** evaluate error when aligned to neighborhood image (every two neighbor images to be aligned)
 * align function format: align(imgX, imgC, width) [suitable for alignPeak and alignShape]
 * @return average error
*/
void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                        int width, std::ofstream &fout);

/** evaluate error when aligned to neighborhood image (every two neighbor images to be aligned)
 * align function format: align(imgX, imgC) [suitable for all align methods]
 * @return error for every image
*/
void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &),
                            const std::string &outPath);

/** evaluate error when aligned to neighborhood image (every two neighbor images to be aligned)
 * align function format: align(imgX, imgC, width) [suitable for alignPeak and alignShape]
 * @return error for every image
*/
void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                            int width, const std::string &outPath);

/** evaluate error which compare similarity between aligned image stack and generated reference image
 * @param dataPath: path for pre-aligned image stack by some method
 * @param refPath: path fot pre-generated reference image by some method
 * @return average error
*/
void evaluateRefImage(const std::string &dataPath, const std::string &refPath, std::ofstream &fout);

#endif //EM_ALIGNEVALUATION_H