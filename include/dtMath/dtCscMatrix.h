/*!
\file       dtCscMatrix.h
\brief      dtMath, Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX_H_
#define DTMATH_DTCSC_MATRIX_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_row, typename m_type> class dtVector;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class dtCscMatrix
{
private:
    int m_elemNum = 0;
    m_type m_elem[m_row * m_col];
    int m_rowIdx[m_row * m_col]; // row indices
    int m_colPtr[m_col + 1];     // column index pointer
    dtCscMatrix(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr);

public:
    dtCscMatrix();
    dtCscMatrix(const dtMatrix<m_row, m_col, m_type> &m);
    dtCscMatrix(const dtCscMatrix &m);
    ~dtCscMatrix() {}

    const m_type *const GetDataAddr() const;
    const int *const GetRowIdx() const;
    const int *const GetColPtr() const;
    uint16_t GetRowSize() const { return m_row; } // size of row
    uint16_t GetColSize() const { return m_col; } // size of colum
    dtVector<m_col, m_type> GetRowVec(const uint16_t idxRow) const;
    dtVector<m_row, m_type> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, dtVector<m_col, m_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, dtVector<m_row, m_type> &v) const;
    dtMatrix<m_row, m_col, m_type> GetDenseMat() const;
    dtCscMatrix<m_col, m_row, m_type> Transpose() const;

    m_type GetNorm() const;
    m_type GetSqNorm() const;

    /* Assignment operators */
    dtCscMatrix &operator=(const dtCscMatrix &m); // matrix  = matrix
    dtCscMatrix &operator*=(const m_type s);      // matrix *= scalar
    dtCscMatrix &operator/=(const m_type s);      // matrix /= scalar

    /* Arithmetic operators */
    dtCscMatrix operator*(const m_type s) const; // matrix * scalar
    dtCscMatrix operator/(const m_type s) const; // matrix / scalar

    dtVector<m_row, m_type> operator*(const dtVector<m_col, m_type> &v) const; // matrix * vector
    dtVector<m_col, m_type> TposeVec(const dtVector<m_row, m_type> &v) const;  // matrix^T * vector

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class dtCscMatrix;
};

} // namespace dtMath

#include "dtCscMatrix.tpp"

#endif // DTMATH_DTCSC_MATRIX_H_
