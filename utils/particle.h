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

#ifndef ALIGNMENT_PARTICLE_H
#define ALIGNMENT_PARTICLE_H


#include "array1D.h"
#include "image2D.h"
#include "noise.h"


template<typename T>
arrayReal<int> coordParticle(const imageReal<T> &img, int size = 80) {
    imageReal<float> imgB = gaussianBlur(img);
    double var = 0;
    int x = 0, y = 0;
    if (size >= img.shape[0] || size >= img.shape[1]) {
        arrayReal<int> coord = {0, 0, imgB.shape[0], imgB.shape[1]};
        return coord;
    }
    for (int i = 0; i < imgB.shape[0] - size; ++i)
        for (int j = 0; j < imgB.shape[1] - size; ++j) {
            imageReal<T> block = imgB.get(i, j, size, size);
            double varT = block.var() / block.mean();
            if (var < varT) {
                var = varT;
                x = i;
                y = j;
            }
        }
    arrayReal<int> coord = {x, y, x + size, y + size};
    return coord;
}

template<typename T>
imageReal<T> getParticle(const imageReal<T> &img, int size = 80) {
    arrayReal<int> coord = coordParticle(img, size);
    imageReal<T> r = img.get(coord[0], coord[1], size, size);
    return r;
}

template<typename T>
imageReal<T> centerParticle(const imageReal<T> &img, int size = 80) {
    arrayReal<int> coord = coordParticle(img, size);
    int dx = coord[0] + size / 2 - img.shape[0] / 2;
    int dy = coord[2] + size / 2 - img.shape[1] / 2;
    imageReal<T> r = shift(img, -dx, -dy);
    return r;
}


#endif //ALIGNMENT_PARTICLE_H
