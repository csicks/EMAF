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

#ifndef ALIGNMENT_MRCHEADER_H
#define ALIGNMENT_MRCHEADER_H

/** reference website:
 * https://bio3d.colorado.edu/imod/doc/mrc_format.txt
 * https://www.ccpem.ac.uk/mrc_format/mrc2014.php */

#include <fstream>
#include <cstring>

#define MODE_INT_8 0
#define MODE_INT_16 1
#define MODE_REAL_32 2
#define MODE_COMPLEX_16 3
#define MODE_COMPLEX_32 4
#define MODE_UINT_16 6

class mrcHeader
{
private:
    int totalSize = 1024;

public:
    int nx;
    int ny;
    int nz;
    int mode;
    int nxStart;
    int nyStart;
    int nzStart;
    int mx;
    int my;
    int mz;
    int cellA[3]{};
    int cellB[3]{};
    int mapC;
    int mapR;
    int mapS;
    float dMin;
    float dMax;
    float dMean;
    int ispg;
    int nSymBT;
    int extra[25]{};
    char extTyp[4];
    int nVersion;
    int origin[3]{};
    char map[4];
    unsigned char machSt[4];
    int rms;
    int nLabl;
    int label[200]{};

public:
    mrcHeader();

    mrcHeader(const mrcHeader &header);

    void read(std::ifstream &mrc);

    void write(std::ofstream &mrc) const;
};

#endif // ALIGNMENT_MRCHEADER_H