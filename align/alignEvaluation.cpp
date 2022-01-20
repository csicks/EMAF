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

#include "alignEvaluation.h"

stackReal<double> readData(const std::string &path) {
    mrcFile mrc = mrcFile();
    mrc.read(path);
    stackReal<double> data = mrc.data;
    return data;
}

std::vector<std::vector<double>> readParams(const std::string &path) {
    std::ifstream in;
    in.open(path, std::ios::in);
    if (!in)
        throw baseException("Error: Unable to open txt files!");
    std::string line;
    std::vector<std::vector<double>> r;
    while (getline(in, line)) {
        std::istringstream sin(line);
        std::vector<double> fields;
        std::string field;
        int count = 0;
        while (getline(sin, field, ' ')) {
            if (count > 0) {
                double number;
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
                       phi (*func)(const imageReal<double> &, const imageReal<double> &), std::ofstream &fout) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    imageReal<double> ref = readData(refPath).pieceGet(0);
    double s = 0;

    for (int i = 0; i < data.shape[0]; ++i) {
        imageReal<double> temp = data.pieceGet(i);
        std::vector<double> params = txt[i];
        phi p = func(temp, ref);
        double angleError = std::abs(params[0] - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        s += angleError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y);
    }

    s /= data.shape[0] * 3;

    std::cout << "RefParams error: " << s << std::endl;
    fout << "RefParams error: " << s << std::endl;
}

void evaluateRefParams(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                       phi (*func)(const imageReal<double> &, const imageReal<double> &, int), int width,
                       std::ofstream &fout) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    imageReal<double> ref = readData(refPath).pieceGet(0);
    double s = 0;

    for (int i = 0; i < data.shape[0]; ++i) {
        imageReal<double> temp = data.pieceGet(i);
        std::vector<double> params = txt[i];
        phi p = func(temp, ref, width);
        double angleError = std::abs(params[0] - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        s += angleError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y);
    }

    s /= data.shape[0] * 3;

    std::cout << "RefParams error: " << s << std::endl;
    fout << "RefParams error: " << s << std::endl;
}

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi (*func)(const imageReal<double> &, const imageReal<double> &),
                           const std::string &outPath) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    imageReal<double> ref = readData(refPath).pieceGet(0);
    std::ofstream fout(outPath, std::ios::out);

    for (int i = 0; i < data.shape[0]; ++i) {
        imageReal<double> temp = data.pieceGet(i);
        std::vector<double> params = txt[i];
        phi p = func(temp, ref);
        double angleError = std::abs(params[0] - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        fout << (angleError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y)) / 3 << std::endl;
    }

    fout.close();
}

void evaluateRefParamsHist(const std::string &dataPath, const std::string &txtPath, const std::string &refPath,
                           phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                           int width, const std::string &outPath) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    imageReal<double> ref = readData(refPath).pieceGet(0);
    std::ofstream fout(outPath, std::ios::out);

    for (int i = 0; i < data.shape[0]; ++i) {
        imageReal<double> temp = data.pieceGet(i);
        std::vector<double> params = txt[i];
        phi p = func(temp, ref, width);
        double angleError = std::abs(params[0] - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        fout << (angleError + std::abs(params[1] - p.x) + std::abs(params[2] - p.y)) / 3 << std::endl;
    }

    fout.close();
}

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi (*func)(const imageReal<double> &, const imageReal<double> &), std::ofstream &fout) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    double s = 0;
    int count = 0;

    for (int i = 0; i + 1 < data.shape[0]; i += 2) {
        imageReal<double> temp1 = data.pieceGet(i);
        imageReal<double> temp2 = data.pieceGet(i + 1);
        std::vector<double> params1 = txt[i];
        std::vector<double> params2 = txt[i + 1];
        phi p = func(temp2, temp1);
        double angleT = params2[0] - params1[0];
        double angleR = angleT / 360 * 2 * M_PI;
        double xT = -params1[2] * std::sin(angleR) - params1[1] * std::cos(angleR) + params2[1];
        double yT = -params1[2] * std::cos(angleR) + params1[1] * std::sin(angleR) + params2[2];
        double angleError = std::abs(angleT - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        s += angleError + std::abs(xT - p.x) + std::abs(yT - p.y);
        ++count;
    }

    s /= count * 3;
    std::cout << "PairParams error: " << s << std::endl;
    fout << "PairParams error: " << s << std::endl;
}

void evaluatePairParams(const std::string &dataPath, const std::string &txtPath,
                        phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                        int width, std::ofstream &fout) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    double s = 0;
    int count = 0;

    for (int i = 0; i + 1 < data.shape[0]; i += 2) {
        imageReal<double> temp1 = data.pieceGet(i);
        imageReal<double> temp2 = data.pieceGet(i + 1);
        std::vector<double> params1 = txt[i];
        std::vector<double> params2 = txt[i + 1];
        phi p = func(temp2, temp1, width);
        double angleT = params2[0] - params1[0];
        double angleR = angleT / 360 * 2 * M_PI;
        double xT = -params1[2] * std::sin(angleR) - params1[1] * std::cos(angleR) + params2[1];
        double yT = -params1[2] * std::cos(angleR) + params1[1] * std::sin(angleR) + params2[2];
        double angleError = std::abs(angleT - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        s += angleError + std::abs(xT - p.x) + std::abs(yT - p.y);
        ++count;
    }

    s /= count * 3;
    std::cout << "PairParams error: " << s << std::endl;
    fout << "PairParams error: " << s << std::endl;
}

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &),
                            const std::string &outPath) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    std::ofstream fout(outPath, std::ios::out);

    for (int i = 0; i + 1 < data.shape[0]; i += 2) {
        imageReal<double> temp1 = data.pieceGet(i);
        imageReal<double> temp2 = data.pieceGet(i + 1);
        std::vector<double> params1 = txt[i];
        std::vector<double> params2 = txt[i + 1];
        phi p = func(temp2, temp1);
        double angleT = params2[0] - params1[0];
        double angleR = angleT / 360 * 2 * M_PI;
        double xT = -params1[2] * std::sin(angleR) - params1[1] * std::cos(angleR) + params2[1];
        double yT = -params1[2] * std::cos(angleR) + params1[1] * std::sin(angleR) + params2[2];
        double angleError = std::abs(angleT - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        fout << (angleError + std::abs(xT - p.x) + std::abs(yT - p.y)) / 3 << std::endl;
    }

    fout.close();
}

void evaluatePairParamsHist(const std::string &dataPath, const std::string &txtPath,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &, int),
                            int width, const std::string &outPath) {
    stackReal<double> data = readData(dataPath);
    std::vector<std::vector<double>> txt = readParams(txtPath);
    std::ofstream fout(outPath, std::ios::out);

    for (int i = 0; i + 1 < data.shape[0]; i += 2) {
        imageReal<double> temp1 = data.pieceGet(i);
        imageReal<double> temp2 = data.pieceGet(i + 1);
        std::vector<double> params1 = txt[i];
        std::vector<double> params2 = txt[i + 1];
        phi p = func(temp2, temp1, width);
        double angleT = params2[0] - params1[0];
        double angleR = angleT / 360 * 2 * M_PI;
        double xT = -params1[2] * std::sin(angleR) - params1[1] * std::cos(angleR) + params2[1];
        double yT = -params1[2] * std::cos(angleR) + params1[1] * std::sin(angleR) + params2[2];
        double angleError = std::abs(angleT - p.angle);
        angleError = std::min(angleError, std::abs(angleError - 360));
        angleError = std::min(angleError, std::abs(angleError - 180));
        fout << (angleError + std::abs(xT - p.x) + std::abs(yT - p.y)) / 3 << std::endl;
    }

    fout.close();
}

void evaluateRefImage(const std::string &dataPath, const std::string &refPath, std::ofstream &fout) {
    stackReal<double> data = readData(dataPath);
    double s = 0;

    imageReal<double> ref = readData(refPath).pieceGet(0);

    for (int i = 0; i < data.shape[0] / 2; ++i) {
        imageReal<double> temp = data.pieceGet(i);
        double value = correntropy(temp, ref);
        s += value;
    }

    s /= data.shape[0];

    std::cout << "RefImage error: " << s << std::endl;
    fout << "RefImage error: " << s << std::endl;
}