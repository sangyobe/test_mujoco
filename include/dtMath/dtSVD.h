/*!
\file       dtSVD.h
\brief      dtMath, Singular Value Decomposition solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD_H_
#define DTMATH_DTSVD_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

/*
This procedure programmed here is based on the method of Golub and Reinsch as given on
pages 134 - 151 of the "Handbook for Automatic Computation vol II - Linear Algebra"
edited by Wilkinson and Reinsch and published by Springer - Verlag, 1971.
*/

/*
This routine decomposes an mxn matrix A, with m >= n,
into a product of the three matrices U, S, and VT,

i.e. A = U S VT,

where U is an mxn matrix whose columns are mutually orthogonal,
      S is an nxn diagonal matrix,
      V is an nxn orthogonal matrix.
      VT denotes the transpose of V.

If m < n, then the procedure may be used for the matrix AT.
The singular values of A are the diagonal elements of the diagonal matrix S
and correspond to the positive square roots of the eigenvalues of the matrix ATA.
*/

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class dtSVD
{
private:
    const int MAX_ITERATION_COUNT = 30;
    m_type m_A[m_row * m_col] = {
        0,
    };
    m_type m_U[(m_row >= m_col) ? m_row * m_col : m_col * m_col] = {
        0,
    };
    m_type m_S[(m_row >= m_col) ? m_col : m_row] = {
        0,
    };
    m_type m_V[(m_row >= m_col) ? m_col * m_col : m_row * m_col] = {
        0,
    };
    m_type m_superDiagonal[(m_row >= m_col) ? m_col : m_row] = {
        0,
    }; // This array is used to store the super-diagonal elements resulting from the Householder reduction of the matrix A to bidiagonal form.
    m_type m_inv[m_row * m_col] = {
        0,
    };
    int8_t m_isOk = 0;

public:
    dtSVD();
    dtSVD(const m_type *element, const size_t n_byte);
    dtSVD(const dtMatrix<m_row, m_col, m_type> &m);
    dtSVD(const dtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute(const m_type *element, const size_t n_byte); // Compute Singular Value Decomposition
    int8_t Compute(const dtMatrix<m_row, m_col, m_type> &m);    // Compute Singular Value Decomposition
    int8_t Compute(const dtMatrix3<m_type, m_row, m_col> &m);   // Compute Singular Value Decomposition
    int8_t IsOk() { return m_isOk; }

    dtMatrix<m_row, m_row, m_type> GetMatrixU() const; // return U matrix, m x m, left singular vectors of A
    dtMatrix<m_row, m_col, m_type> GetMatrixS() const; // return S matrix, m x n, singular values of A
    dtMatrix<m_col, m_col, m_type> GetMatrixV() const; // return V matrix, n x n, right singular vectors of A

    template <uint16_t col>
    int8_t Solve(const dtMatrix<m_row, col, m_type> &b, dtMatrix<m_col, col, m_type> &x, m_type tolerance = 0); // Solve x = (USVT)^-1 * b
    int8_t Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x, m_type tolerance = 0);           // Solve x = (USVT)^-1 * b
    template <uint16_t col>
    dtMatrix<m_col, col, m_type> Solve(const dtMatrix<m_row, col, m_type> &b, int8_t *isOk = nullptr, m_type tolerance = 0); // Solve x = (USVT)^-1 * b
    dtVector<m_col, m_type> Solve(const dtVector<m_row, m_type> &b, int8_t *isOk = nullptr, m_type tolerance = 0);           // Solve x = (USVT)^-1 * b

    int8_t Inverse(dtMatrix<m_col, m_row, m_type> &inv, m_type tolerance = 0);            // Inverse matrix of USVT matrix
    int8_t Inverse(dtMatrix3<m_type, m_col, m_row> &inv, m_type tolerance = 0);           // Inverse matrix of USVT matrix
    dtMatrix<m_col, m_row, m_type> Inverse(int8_t *isOk = nullptr, m_type tolerance = 0); // Inverse matrix of USVT matrix

    int8_t InverseArray(m_type *inv, m_type tolerance = 0);             // Inverse array of USVT matrix
    m_type *InverseArray(int8_t *isOk = nullptr, m_type tolerance = 0); // Inverse array of USVT matrix

private:
    void HouseholdersReductionToBidiagonalForm_Mxn();
    int8_t GivensReductionToDiagonalForm_Mxn();
    void SortByDecreasingSingularValues_Mxn();

    void HouseholdersReductionToBidiagonalForm_mxN();
    int8_t GivensReductionToDiagonalForm_mxN();
    void SortByDecreasingSingularValues_mxN();
};

} // namespace dtMath

#include "dtSVD.tpp"

#endif // DTMATH_DTSVD_H_
