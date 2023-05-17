/*!
\file       dtPartialPivLU.h
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTPARTIAL_PIV_LU_H_
#define DTMATH_DTPARTIAL_PIV_LU_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class dtPartialPivLU
{
private:
    m_type m_elem[m_row * m_col];
    m_type m_inv[m_row * m_col];
    int m_pivot[m_row];
    int8_t m_isOk;

public:
    dtPartialPivLU();
    dtPartialPivLU(const m_type *element, const size_t n_byte);
    dtPartialPivLU(const dtMatrix<m_row, m_col, m_type> &m);
    dtPartialPivLU(const dtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute();                                           // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const m_type *element, const size_t n_byte); // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const dtMatrix<m_row, m_col, m_type> &m);    // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const dtMatrix3<m_type, m_row, m_col> &m);   // Compute Lower/Upper Triangular Matrix, Doolittle form
    m_type Determinant();
    int8_t IsOk() { return m_isOk; }

    dtMatrix<m_row, m_col, m_type> GetMatrix() const;  // return matrix A including L/U matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixL() const; // return Lower Triangular matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixU() const; // return Upper Triangular matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixP() const; // return Permutation matrix

    int8_t Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);              // Solve x = (LU)^-1 * b
    dtVector<m_col, m_type> Solve(const dtVector<m_row, m_type> &b, int8_t *isOk = nullptr); // Solve x = (LU)^-1 * b

    int8_t Inverse(dtMatrix<m_row, m_col, m_type> &inv);            // Inverse matrix of LU matrix
    int8_t Inverse(dtMatrix3<m_type, m_row, m_col> &inv);           // Inverse matrix of LU matrix
    dtMatrix<m_row, m_col, m_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LU matrix

    int8_t InverseArray(m_type *inv);             // Inverse array of LU matrix
    m_type *InverseArray(int8_t *isOk = nullptr); // Inverse array of LU matrix
};

} // namespace dtMath

#include "dtPartialPivLU.tpp"

#endif // DTMATH_DTPARTIAL_PIV_LU_H_
