/*!
\file       dtLDLT.h
\brief      dtMath, Cholesky decomposition(L*D*L^T form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLDLT_H_
#define DTMATH_DTLDLT_H_

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
class dtLDLT
{
private:
    m_type m_elem[m_row * m_col];
    int8_t m_isOk;

public:
    dtLDLT();
    dtLDLT(const m_type *element, const size_t n_byte);
    dtLDLT(const dtMatrix<m_row, m_col, m_type> &m);
    dtLDLT(const dtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute();                                           // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const m_type *element, const size_t n_byte); // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const dtMatrix<m_row, m_col, m_type> &m);    // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const dtMatrix3<m_type, m_row, m_col> &m);   // Compute Cholesky Decomposition, L*D*L^T form
    int8_t IsOk() { return m_isOk; }

    dtMatrix<m_row, m_col, m_type> GetMatrix() const;  // return matrix A including L/D/U matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixL() const; // return Lower Triangular matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixD() const; // return Diagonal matrix
    dtMatrix<m_row, m_col, m_type> GetMatrixU() const; // return Upper Triangular matrix

    template <uint16_t col>
    int8_t Solve(const dtMatrix<m_row, col, m_type> &b, dtMatrix<m_col, col, m_type> &x); // Solve x = (LDU)^-1 * b
    int8_t Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);           // Solve x = (LDU)^-1 * b
    template <uint16_t col>
    dtMatrix<m_col, col, m_type> Solve(const dtMatrix<m_row, col, m_type> &b, int8_t *isOk = nullptr); // Solve x = (LDU)^-1 * b
    dtVector<m_col, m_type> Solve(const dtVector<m_row, m_type> &b, int8_t *isOk = nullptr);           // Solve x = (LDU)^-1 * b

    int8_t Inverse(dtMatrix<m_row, m_col, m_type> &inv);            // Inverse matrix of LDU matrix
    int8_t Inverse(dtMatrix3<m_type, m_row, m_col> &inv);           // Inverse matrix of LDU matrix
    dtMatrix<m_row, m_col, m_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LDU matrix

    int8_t InverseArray(m_type *inv);             // Inverse array of LDU matrix
    m_type *InverseArray(int8_t *isOk = nullptr); // Inverse array of LDU matrix
};

} // namespace dtMath

#include "dtLDLT.tpp"

#endif // DTMATH_DTLDLT_H_
