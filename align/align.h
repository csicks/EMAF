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

#ifndef EM_ALIGN_H
#define EM_ALIGN_H

/** @file
 * this file mainly contains functions relevant to 2D image alignment which includes pair-wise alignment, alignment
 * for image stack and some support functions
*/

#include "image2D.h"
#include "fft.h"
#include "stack3D.h"
#include "polar.h"
#include "fileStack.h"

#define ALIGN_TH 1 // threshold for alignment refinement, which is used to decide whether to stop the iteration

/// structure containing one rotation angle and two translations
struct phi {
    double angle = 0;
    double x = 0;
    double y = 0;

    phi operator+(phi p) const {
        phi r;
        r.angle = angle + p.angle;
        r.x = x + p.x;
        r.y = y + p.y;
        return r;
    }

    phi() = default;

    phi(const std::initializer_list<double> &initList) {
        const double *p = initList.begin();
        angle = *p++;
        x = *p++;
        y = *p;
    }
};

/** common used variable explanation:
 * @param imgX: input image for align function
 * @param imgC: reference image for align function
 * @param stk: input image/file stack
 * @param ref: given reference image for align refinement
 * @param func: given align method [function pointer]
 * @param updateFlag: whether to update reference image during the refinement
 * @param fineFlag: whether to apply fine search for transformation parameters
 * @param info: some head information contained in XMD file
 * @param xmdName: path for output XMD file if needed
*/

/// only align for translations
phi shiftFFT(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// only align for rotation
phi rotateFFT(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// align for translation used in XMIPP
phi shiftFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// align for rotation used in XMIPP
phi rotateFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// convert image to Fourier space and then convert to polar coordinates
imageComplex polarFFT(const imageReal<double> &img);

/** align for both rotation and translation
 * @version 1
*/
phi alignFFT(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// given one rough rotation angle, try to consider symmetry and give out the best rotation angle
phi angleAll(double angle, const imageReal<double> &imgX, const imageReal<double> &imgC, bool fineFlag = false);

/** align for both rotation and translation
 * @version 2
*/
phi alignFFT2(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// align for both rotation and translation, transversing all translations
phi alignShiftPolar(const imageReal<double> &imgX, const imageReal<double> &imgC, int maxShiftX = 4, int maxShiftY = 4);

/// align for both rotation and translation, used in XMIPP
phi alignFFTX(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// apply image transformation according to given transformation parameters
imageReal<double> applyPhi(const imageReal<double> &img, phi p);

/// apply inverse image transformation according to given transformation parameters
imageReal<double> inversePhi(const imageReal<double> &img, phi p);

/** generate reference image for given image stack, SPIDER method
 * align function format: align(imgX, imgC) [suitable for all align methods]
*/
imageReal<double>
reference(const stackReal<double> &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &));

/** generate reference image for given image stack, SPIDER method
 * align function format: align(imgX, imgC, width, fineFlag) [suitable for alignPeak and alignShape]
*/
imageReal<double>
reference(const stackReal<double> &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool),
          int width, bool fineFlag = false);

/** generate reference image for given image stack, SPIDER method
 * align function format: align(imgX, imgC, outWidth, inWidth, th, fineFlag) [suitable for alignShape]
*/
imageReal<double> reference(const stackReal<double> &stk,
                            phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double, bool),
                            int outWidth, int inWidth = -1, double th = -1, bool fineFlag = false);

/** generate reference image for file stack, SPIDER method
 * align function format: align(imgX, imgC) [suitable for all align methods]
*/
imageReal<double> reference(const fileStack &stk, phi (*func)(const imageReal<double> &, const imageReal<double> &));

/// try to move particle to the center of image
imageReal<double> centerImage(const imageReal<double> &img);

/// generate reference image for image stack, XMIPP method
imageReal<double> initReference(const stackReal<double> &stk);

/// generate reference image for file stack, XMIPP method
imageReal<double> initReference(const fileStack &stk);

/** apply alignment to image stack
 * align function format: align(imgX, imgC) [suitable for all align methods]
*/
void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &), bool updateFlag = true);

/** apply alignment to image stack
 * align function format: align(imgX, imgC, width, fineFlag) [suitable for alignPeak and alignShape]
*/
void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool), int width,
                bool fineFlag = false, bool updateFlag = true);

/** apply alignment to image stack
 * align function format: align(imgX, imgC, outWidth, inWidth, th, fineFlag) [suitable for alignShape]
*/
void refinement(stackReal<double> &stk, imageReal<double> &ref,
                phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double, bool), int outWidth,
                int inWidth = -1, double th = -1, bool fineFlag = false, bool updateFlag = true);

/** apply alignment to image stack, with XMD file output
 * align function format: align(imgX, imgC) [suitable for all align methods]
*/
void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &), std::vector<std::string> info,
             const std::string &xmdName);

/** apply alignment to image stack, with XMD file output
 * align function format: align(imgX, imgC, width, fineFlag) [suitable for alignPeak and alignShape]
*/
void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &, int, bool),
             std::vector<std::string> info, const std::string &xmdName, int width, bool fineFlag = false);

/** apply alignment to image stack, with XMD file output
 * align function format: align(imgX, imgC, outWidth, inWidth, th, fineFlag) [suitable for alignShape]
*/
void outInfo(stackReal<double> &stk, imageReal<double> &ref,
             phi (*func)(const imageReal<double> &, const imageReal<double> &, int, int, double, bool),
             std::vector<std::string> info, const std::string &xmdName, int outWidth, int inWidth = -1, double th = -1,
             bool fineFlag = false);

/** apply alignment to file stack
 * align function format: align(imgX, imgC) [suitable for all align methods]
*/
void
refinement(fileStack &stk, imageReal<double> &ref, phi (*func)(const imageReal<double> &, const imageReal<double> &));

/** apply alignment to image stack, align function specified in definition
 * this function have two stks with one involves in alignment but the other only applys the transformation
 * @deprecated
*/
void refinementWithAux(stackReal<double> &stk, stackReal<double> &stkAux, imageReal<double> &ref);

/** a pair-wise align method in a paper (may take some time to find out)
 * @deprecated
*/
phi alignPower(const imageReal<double> &imgX, const imageReal<double> &imgC);

/** EMAF-LM align method, which use local maximum in Fourier space as feature
 * parameters specified in definition
*/
phi alignPeak(const imageReal<double> &imgX, const imageReal<double> &imgC);

/// estimate high-pass filter width for alignPeak method
int circleEstimatePeak(const stackReal<double> &stk, bool fineFlag = false);

/// estimate low-pass filter width for alignShape method
int circleEstimateShape(const stackReal<double> &stk, int inWidth = -1, double th = -1, bool fineFlag = false);

/** EMAF-LM align method, which use local maximum in Fourier space as feature
 * @param width: width for high-pass filter
*/
phi alignPeak(const imageReal<double> &imgX, const imageReal<double> &imgC, int width, bool fineFlag = false);

/** EMAF-MD align method, which use moment direction in Fourier space as feature
 * parameters specified in definition
 * */
phi alignShape(const imageReal<double> &imgX, const imageReal<double> &imgC);

/** EMAF-MD align method, which use moment direction in Fourier space as feature
 * @param outWidth: width for low-pass filter
 * @param inWidth: width for high-pass filter
 * @param th: threshold value for threshold cutting
*/
phi
alignShape(const imageReal<double> &imgX, const imageReal<double> &imgC, int outWidth, int inWidth = -1, double th = -1,
           bool fineFlag = false);

/// pair-wise alignment method, which find and extract particle firstly and then apply alignment
phi alignParticle(const imageReal<double> &imgX, const imageReal<double> &imgC);

#endif //EM_ALIGN_H
