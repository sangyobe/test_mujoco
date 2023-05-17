/*!
\file       dtMatrix3.h
\brief      dtMath, 3x3 Matrix class, lighter and faster than general matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX3_H_
#define DTMATH_DTMATRIX3_H_

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

template <uint16_t m_size, typename m_type> class dtCommaInit;
template <uint16_t m_row, typename m_type> class dtVector;
template <typename m_type, uint16_t m_row> class dtVector3;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtNoPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtPartialPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtLLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtLDLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtQR;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtSVD;
template <typename type, uint16_t row, uint16_t col> class dtRotation;

template <typename m_type = float, uint16_t m_row = 3, uint16_t m_col = 3>
class dtMatrix3
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    dtMatrix3(const m_type *element);

public:
    dtMatrix3();
    dtMatrix3(const m_type *element, const size_t n_byte);
    dtMatrix3(const char c, const m_type *element, const size_t n_byte);
    dtMatrix3(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    dtMatrix3(const dtMatrix3 &m);
    dtMatrix3(const dtRotation<m_type, m_row, m_col> &m);
    dtMatrix3(const dtMatrix<m_row, m_col, m_type> &m);
    ~dtMatrix3() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const m_type d1, const m_type d2, const m_type d3);
    void SetDiagonal(const m_type *element, const size_t n_byte);
    void SetDiagonal(const dtVector<m_row, m_type> &v);
    void SetDiagonal(const dtVector3<m_type, m_row> &v);
    void SetFill(const m_type value);
    void SetElement(const m_type *element, const size_t n_byte);
    void SetElement(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    void SetElement(const dtMatrix3 &m);
    void SetElement(const dtRotation<m_type, m_row, m_col> &m);
    void SetElement(const dtMatrix<m_row, m_col, m_type> &m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const dtVector<col, m_type> &v);
    void SetRowVec(const uint16_t idxRow, const dtVector3<m_type, m_col> &v);
    void SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const dtVector<row, m_type> &v);
    void SetColVec(const uint16_t idxCol, const dtVector3<m_type, m_row> &v);
    void SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type *const GetElementsAddr() const;
    dtVector3<m_type, m_col> GetRowVec(const uint16_t idxRow) const;
    dtVector3<m_type, m_row> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, dtVector3<m_type, m_col> &v) const;
    int8_t GetColVec(const uint16_t idxCol, dtVector3<m_type, m_row> &v) const;
    dtMatrix3 Transpose() const;

    m_type Trace() const;
    m_type GetNorm() const;   // Frobenius Norm (Euclidean norm)
    m_type GetSqNorm() const; // Squared Frobenius Norm (Euclidean norm)

    dtNoPivLU<m_row, m_col, m_type> NoPivLU() const;
    dtPartialPivLU<m_row, m_col, m_type> PartialPivLU() const;
    dtLLT<m_row, m_col, m_type> LLT() const;
    dtLDLT<m_row, m_col, m_type> LDLT() const;
    dtQR<m_row, m_col, m_type> QR() const;
    dtSVD<m_row, m_col, m_type> SVD() const;

    dtMatrix3 Inv(int8_t *isOk = nullptr) const;
    dtMatrix3 PInv(int8_t *isOk = nullptr, m_type tolerance = 0) const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    dtMatrix3 &operator=(const dtMatrix3 &m);                         // matrix  = matrix
    dtMatrix3 &operator+=(const dtMatrix3 &m);                        // matrix += matrix
    dtMatrix3 &operator-=(const dtMatrix3 &m);                        // matrix -= matrix
    dtMatrix3 &operator=(const dtRotation<m_type, m_row, m_col> &m);  // matrix  = matrix
    dtMatrix3 &operator+=(const dtRotation<m_type, m_row, m_col> &m); // matrix += matrix
    dtMatrix3 &operator-=(const dtRotation<m_type, m_row, m_col> &m); // matrix -= matrix
    dtMatrix3 &operator=(const dtMatrix<m_row, m_col, m_type> &m);    // matrix  = matrix
    dtMatrix3 &operator+=(const dtMatrix<m_row, m_col, m_type> &m);   // matrix += matrix
    dtMatrix3 &operator-=(const dtMatrix<m_row, m_col, m_type> &m);   // matrix -= matrix
    dtMatrix3 &operator=(const m_type s);                             // matrix  = scalar, all elements set scalar
    dtMatrix3 &operator+=(const m_type s);                            // matrix += scalar, matrix(i) += scalar
    dtMatrix3 &operator-=(const m_type s);                            // matrix -= scalar, matrix(i) -= scalar
    dtMatrix3 &operator*=(const m_type s);                            // matrix *= scalar
    dtMatrix3 &operator/=(const m_type s);                            // matrix /= scalar
    dtCommaInit<m_row * m_col, m_type> operator<<(const m_type s);    // Init first matrix elements

    /* Arithmetic operators */
    dtMatrix3 operator-() const;                                          // minus sign
    dtMatrix3 operator+(const dtMatrix3 &m) const;                        // matrix + matrix
    dtMatrix3 operator-(const dtMatrix3 &m) const;                        // matrix - matrix
    dtMatrix3 operator+(const dtRotation<m_type, m_row, m_col> &m) const; // matrix + matrix
    dtMatrix3 operator-(const dtRotation<m_type, m_row, m_col> &m) const; // matrix - matrix
    dtMatrix3 operator+(const dtMatrix<m_row, m_col, m_type> &m) const;   // matrix + matrix
    dtMatrix3 operator-(const dtMatrix<m_row, m_col, m_type> &m) const;   // matrix - matrix
    dtMatrix3 operator+(const m_type s) const;                            // matrix + scalar, matrix(i) + scalar
    dtMatrix3 operator-(const m_type s) const;                            // matrix - scalar, matrix(i) - scalar
    dtMatrix3 operator*(const m_type s) const;                            // matrix * scalar
    dtMatrix3 operator/(const m_type s) const;                            // matrix / scalar

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<m_col, col, m_type> &m) const; // matrix * matrix
    dtMatrix3 operator*(const dtMatrix3 &m) const;                                       // matrix * matrix
    dtMatrix3 operator*(const dtRotation<m_type, m_row, m_col> &m) const;                // matrix * rotation
    dtVector<m_row, m_type> operator*(const dtVector<m_col, m_type> &v) const;           // matrix * vector
    dtVector3<m_type, m_row> operator*(const dtVector3<m_type, m_col> &v) const;         // matrix * vector
    dtMatrix3 operator&(const dtVector<m_col, m_type> &v) const;                         // matrix * [v]x, []x is skew-symmetric matrix
    dtMatrix3 operator&(const dtVector3<m_type, m_col> &v) const;                        // matrix * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const dtMatrix3 &m) const;                        // (true or false) matrix == matrix
    bool operator!=(const dtMatrix3 &m) const;                        // (true or false) matrix != matrix
    bool operator==(const dtRotation<m_type, m_row, m_col> &m) const; // (true or false) matrix == matrix
    bool operator!=(const dtRotation<m_type, m_row, m_col> &m) const; // (true or false) matrix != matrix
    bool operator==(const dtMatrix<m_row, m_col, m_type> &m) const;   // (true or false) matrix == matrix
    bool operator!=(const dtMatrix<m_row, m_col, m_type> &m) const;   // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class dtMatrix3;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;

    template <uint16_t row, uint16_t col, typename type> friend class dtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtQR;
    template <uint16_t row, uint16_t col, typename type> friend class dtSVD;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend dtMatrix3<type, row, col> operator*(const type s, const dtMatrix3<type, row, col> &m); // scalar * matrix
};

} // namespace dtMath

#include "dtMatrix3.tpp"

#endif // DTMATH_DTMATRIX3_H_
