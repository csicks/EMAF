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

#include "fileStack.h"


std::vector<std::string> split(const std::string &line, char symbol) {
    std::istringstream sin(line);
    std::vector<std::string> fields;
    std::string field;
    while (getline(sin, field, symbol)) {
        fields.push_back(field);
    }
    return fields;
}


fileStack::fileStack(const std::string &aux, const std::string &auxB) {
    if (access(auxB.c_str(), 0) == 0)
        system(("rm -r " + auxB).c_str());
    int mdB = mkdir(auxB.c_str(), S_IRWXU);
    if (mdB == -1)
        throw baseException("Error: Unable to create aux folder for file stacks!");
    auxFolder = aux;
    auxFolderBack = auxB;
}


fileStack::fileStack(const std::string &path, const std::string &extension, const std::string &aux,
                     const std::string &auxB) {
    if (access(aux.c_str(), 0) == 0)
        system(("rm -r " + aux).c_str());
    if (access(auxB.c_str(), 0) == 0)
        system(("rm -r " + auxB).c_str());
    int md = mkdir(aux.c_str(), S_IRWXU);
    int mdB = mkdir(auxB.c_str(), S_IRWXU);
    if (md == -1 || mdB == -1)
        throw baseException("Error: Unable to create aux folder for file stacks!");
    auxFolder = aux;
    auxFolderBack = auxB;

    glob_t g;
    std::string pathAll = path + extension;
    const char *pattern = pathAll.c_str();

    if (glob(pattern, GLOB_ERR, nullptr, &g) != 0)
        throw baseException("Error: Failed to load from folder!");

    for (int i = 0; i < g.gl_pathc; i++) {
        mrcFile mrc = mrcFile();
        mrc.read(g.gl_pathv[i]);
        std::vector<std::string> sg = split(g.gl_pathv[i], '/');
        std::string fileName = sg[sg.size() - 1];
        mrc.write(auxFolder + fileName);
        names.emplace_back(auxFolder + fileName);
    }
    globfree(&g);

    mrcFile mrc = mrcFile();
    mrc.read(names[0]);
    shape[0] = names.size();
    shape[1] = mrc.data.shape[1];
    shape[2] = mrc.data.shape[2];
}


float fileStack::get(int i, int j, int k) const {
    imageReal<float> img = pieceGet(i);
    return img.get(j, k);
}


imageReal<float> fileStack::pieceGet(int index) const {
    mrcFile mrc = mrcFile();
    mrc.read(names[index]);
    imageReal<float> img = mrc.data.pieceGet(0);
    return img;
}


void fileStack::pieceSet(int index, const imageReal<float> &img) {
    mrcFile mrc_new = mrcFile();
    stackReal<float> stk = image2Stack(img);
    mrc_new.setData(stk);
    std::string pathOri = names[index];
    std::vector<std::string> sg = split(pathOri, '/');
    std::string fileName = sg[sg.size() - 1];
    std::string path = auxFolder + fileName;
    mrc_new.write(path);
    names[index] = path;
}


void fileStack::append(const std::string &path) {
    if (shape[0] == 0) {
        mrcFile mrc = mrcFile();
        mrc.read(path);
        shape[1] = mrc.data.shape[1];
        shape[2] = mrc.data.shape[2];
    }
    shape[0]++;
    names.push_back(path);
}


void fileStack::append(const imageReal<float> &img) {
    if (shape[0] == 0) {
        shape[1] = img.shape[0];
        shape[2] = img.shape[1];
    }
    shape[0]++;
    mrcFile mrc = mrcFile();
    stackReal<float> stk = image2Stack(img);
    mrc.setData(stk);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0, 1);
    auto gr = static_cast<float>(distribution(generator));
    int randomLabel = static_cast<int>(gr * 10000000);
    randomLabel = randomLabel < 0 ? -randomLabel : randomLabel;
    std::string ss = std::to_string(randomLabel);
    int width = 10;
    while (ss.size() < width) {
        ss = "0" + ss;
    }
    std::string name = auxFolderBack + ss + ".mrcs";
    std::string nameNew = auxFolder + ss + ".mrcs";

    mrc.write(name);
    names.emplace_back(nameNew);
}