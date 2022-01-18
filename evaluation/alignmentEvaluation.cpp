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

#include "alignmentEvaluation.h"

stackReal<float> readData(const std::string &path)
{
    mrcFile mrc = mrcFile();
    mrc.read(path);
    stackReal<float> data = mrc.data;
    return data;
}

std::vector<std::vector<float>> readParams(const std::string &path)
{
    std::ifstream in;
    in.open(path, std::ios::in);
    if (!in)
        throw baseException("Error: Unable to open txt file!");
    std::string line;
    std::vector<std::vector<float>> r;
    while (getline(in, line))
    {
        std::istringstream sin(line);
        std::vector<float> fields;
        std::string field;
        int count = 0;
        while (getline(sin, field, ' '))
        {
            if (count > 0)
            {
                float number;
                std::istringstream stream(field);
                stream >> number;
                fields.push_back(number);
            }
            ++count;
        }
        r.push_back(fields);
    }
    return r;
}

void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    imageReal<float> ref = readData(refPath).pieceGet(0);
    float s = 0;

    for (int i = 0; i < data.shape[0]; ++i)
    {
        imageReal<float> temp = data.pieceGet(i);
        std::vector<float> params = txt[i];
        phi p = func(temp, ref);
        double angError = std::abs(params[0] - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        s += angError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y);
    }

    s /= data.shape[0] * 3;

    std::cout << "RefParams error: " << s << std::endl;
    fout << "RefParams error: " << s << std::endl;
}

void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    imageReal<float> ref = readData(refPath).pieceGet(0);
    float s = 0;

    for (int i = 0; i < data.shape[0]; ++i)
    {
        imageReal<float> temp = data.pieceGet(i);
        std::vector<float> params = txt[i];
        phi p = func(temp, ref, width);
        double angError = std::abs(params[0] - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        s += angError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y);
    }

    s /= data.shape[0] * 3;

    std::cout << "RefParams error: " << s << std::endl;
    fout << "RefParams error: " << s << std::endl;
}

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    imageReal<float> ref = readData(refPath).pieceGet(0);
    std::ofstream outFile(outPath, std::ios::out);

    for (int i = 0; i < data.shape[0]; ++i)
    {
        imageReal<float> temp = data.pieceGet(i);
        std::vector<float> params = txt[i];
        phi p = func(temp, ref);
        double angError = std::abs(params[0] - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        outFile << (angError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y)) / 3 << std::endl;
    }

    outFile.close();
}

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    imageReal<float> ref = readData(refPath).pieceGet(0);
    std::ofstream outFile(outPath, std::ios::out);

    for (int i = 0; i < data.shape[0]; ++i)
    {
        imageReal<float> temp = data.pieceGet(i);
        std::vector<float> params = txt[i];
        phi p = func(temp, ref, width);
        double angError = std::abs(params[0] - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        outFile << (angError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y)) / 3 << std::endl;
    }

    outFile.close();
}

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), std::ofstream &fout)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    float s = 0;
    int count = 0;

    for (int i = 0; i + 1 < data.shape[0]; i += 2)
    {
        imageReal<float> temp1 = data.pieceGet(i);
        imageReal<float> temp2 = data.pieceGet(i + 1);
        std::vector<float> params1 = txt[i];
        std::vector<float> params2 = txt[i + 1];
        phi p = func(temp2, temp1);
        float angT = params2[0] - params1[0];
        float angR = angT / 360 * 2 * M_PI;
        float xT = -params1[2] * std::sin(angR) - params1[1] * std::cos(angR) + params2[1];
        float yT = -params1[2] * std::cos(angR) + params1[1] * std::sin(angR) + params2[2];
        double angError = std::abs(angT - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        s += angError + std::abs(xT - p.x) + std::abs(yT - p.y);
        ++count;
    }

    s /= count * 3;
    std::cout << "PairParams error: " << s << std::endl;
    fout << "PairParams error: " << s << std::endl;
}

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, std::ofstream &fout)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    float s = 0;
    int count = 0;

    for (int i = 0; i + 1 < data.shape[0]; i += 2)
    {
        imageReal<float> temp1 = data.pieceGet(i);
        imageReal<float> temp2 = data.pieceGet(i + 1);
        std::vector<float> params1 = txt[i];
        std::vector<float> params2 = txt[i + 1];
        phi p = func(temp2, temp1, width);
        float angT = params2[0] - params1[0];
        float angR = angT / 360 * 2 * M_PI;
        float xT = -params1[2] * std::sin(angR) - params1[1] * std::cos(angR) + params2[1];
        float yT = -params1[2] * std::cos(angR) + params1[1] * std::sin(angR) + params2[2];
        double angError = std::abs(angT - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        s += angError + std::abs(xT - p.x) + std::abs(yT - p.y);
        ++count;
    }

    s /= count * 3;
    std::cout << "PairParams error: " << s << std::endl;
    fout << "PairParams error: " << s << std::endl;
}

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC), const std::string &outPath)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    std::ofstream outFile(outPath, std::ios::out);

    for (int i = 0; i + 1 < data.shape[0]; i += 2)
    {
        imageReal<float> temp1 = data.pieceGet(i);
        imageReal<float> temp2 = data.pieceGet(i + 1);
        std::vector<float> params1 = txt[i];
        std::vector<float> params2 = txt[i + 1];
        phi p = func(temp2, temp1);
        float angT = params2[0] - params1[0];
        float angR = angT / 360 * 2 * M_PI;
        float xT = -params1[2] * std::sin(angR) - params1[1] * std::cos(angR) + params2[1];
        float yT = -params1[2] * std::cos(angR) + params1[1] * std::sin(angR) + params2[2];
        double angError = std::abs(angT - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        outFile << (angError + std::abs(xT - p.x) + std::abs(yT - p.y)) / 3 << std::endl;
    }

    outFile.close();
}

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi func(const imageReal<float> &imgX, const imageReal<float> &imgC, int width), int width, const std::string &outPath)
{
    stackReal<float> data = readData(dataPath);
    std::vector<std::vector<float>> txt = readParams(txtPath);
    std::ofstream outFile(outPath, std::ios::out);

    for (int i = 0; i + 1 < data.shape[0]; i += 2)
    {
        imageReal<float> temp1 = data.pieceGet(i);
        imageReal<float> temp2 = data.pieceGet(i + 1);
        std::vector<float> params1 = txt[i];
        std::vector<float> params2 = txt[i + 1];
        phi p = func(temp2, temp1, width);
        float angT = params2[0] - params1[0];
        float angR = angT / 360 * 2 * M_PI;
        float xT = -params1[2] * std::sin(angR) - params1[1] * std::cos(angR) + params2[1];
        float yT = -params1[2] * std::cos(angR) + params1[1] * std::sin(angR) + params2[2];
        double angError = std::abs(angT - p.ang);
        angError = std::min(angError, std::abs(angError - 360));
        angError = std::min(angError, std::abs(angError - 180));
        outFile << (angError + std::abs(xT - p.x) + std::abs(yT - p.y)) / 3 << std::endl;
    }

    outFile.close();
}

void evaluateRefImage(const std::string &dataPath, const std::string &refPath, std::ofstream &fout)
{
    stackReal<float> data = readData(dataPath);
    float s = 0;

    imageReal<float> ref = readData(refPath).pieceGet(0);

    for (int i = 0; i < data.shape[0] / 2; ++i)
    {
        imageReal<float> temp = data.pieceGet(i);
        float value = correntropy(temp, ref);
        s += value;
    }

    s /= data.shape[0];

    std::cout << "RefImage error: " << s << std::endl;
    fout << "RefImage error: " << s << std::endl;
}