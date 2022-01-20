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

#ifndef EM_MATRIX_H
#define EM_MATRIX_H

/** @file
 * this file contains data structure for two-dimensional real matrix
 * TODO: may add support to complex case
 * TODO: may replace get/set with () and remove [] (which will be supported in C++23)
*/

#include "array1D.h"

/** @brief
 * class for two-dimensional real matrix, whose size could be changed
*/
class matrix {
public:
    int shape[2]{}; // shape of matrix
    arrayReal<double> data; // data stored in one-dimensional array

public:
    matrix();

    matrix(const matrix &matrix);

    explicit matrix(const int s[2]);

    matrix(const std::initializer_list<double> &initList, const int s[2]);

    matrix(const arrayReal<double> &d, const int *s);

    matrix(const double *d, const int s[2]);

    matrix &operator=(const matrix &matrix);

    /// @return value at position (i,j), i.e. M[i,j]
    double get(int i, int j) const;

    /// set M[i,j] to value
    void set(int i, int j, double value);

    /** @attention assignment value, i.e. matrix[i]=ary, doesn't work!
     * @return the ith row of matrix
    */
    arrayReal<double> operator[](int index) const;

    /// @return the ith row of matrix
    arrayReal<double> row(int index) const;

    /// @return the ith column of matrix
    arrayReal<double> column(int index) const;

    /// remove the ith row of matrix
    matrix removeRow(int index) const;

    /// remove the ith column of matrix
    matrix removeColumn(int index) const;

    matrix operator+(const matrix &matrixIn) const;

    matrix operator-(const matrix &matrixIn) const;

    matrix operator*(const matrix &matrixIn) const;

    template<typename T>
    matrix operator+(T value) const {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data + value;
        return matrixOut;
    }

    template<typename T>
    matrix operator-(T value) const {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data - value;
        return matrixOut;
    }

    template<typename T>
    matrix operator*(T value) const {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data * value;
        return matrixOut;
    }

    template<typename T>
    matrix operator/(T value) const {
        matrix matrixOut = matrix(shape);
        matrixOut.data = data / value;
        return matrixOut;
    }

    /// @return determinant of matrix
    double determinant() const;

    /// @return transpose matrix of this one
    matrix transpose() const;

    /// @return adjoint matrix of this one
    matrix adjoint() const;

    /** @return inverse matrix of this one
     * @implements @code
     *             $M^{-1}=\frac{adjM}{|M|}$
     *             @endcode
     * TODO: try other algorithms to speed up inverse matrix computation
    */
    matrix inverse() const;

    /// print matrix
    void print() const;
};

/// matrix used for two-dimensional image transformation, shape of [3,3]
class matrixTransform : public matrix {

public:
    using matrix::matrix;

    /** generate transformation matrix from rotation angle and two translations
     * @code
     * $$\left[\begin{array}{ccc}
     * \cos(\alpha) & \sin(\alpha) & x \\
     * -\sin(\alpha) & \cos(\alpha) & y \\
     * 0 & 0 & 1
     * \end{array}\right]$$
     * @endcode
    */
    matrixTransform(double angle, int x, int y);

    matrixTransform(const std::initializer_list<double> &initList);

    explicit matrixTransform(const matrix &m);

    matrixTransform &operator=(const matrixTransform &m);
};

#endif //EM_MATRIX_H