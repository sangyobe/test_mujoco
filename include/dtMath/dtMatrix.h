/*!
\file       dtMatrix.h
\brief      dtMath, General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX_H_
#define DTMATH_DTMATRIX_H_

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
template <typename m_type, uint16_t m_row> class dtVector4;
template <typename m_type, uint16_t m_row> class dtVector6;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtRotation;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtTransform;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtNoPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtPartialPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtLLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtLDLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtQR;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtSVD;

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class dtMatrix
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    dtMatrix(const m_type *element);

public:
    dtMatrix();
    dtMatrix(const m_type *element, const size_t n_byte);
    dtMatrix(const char c, const m_type *element, const size_t n_byte = m_row * m_col); // JhJo(230424) : modified to handle default n_byte as matrix size.
    dtMatrix(const dtMatrix &m);
    ~dtMatrix() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const m_type *element, const size_t n_byte);
    void SetFill(const m_type value);
    void SetElement(const m_type *element, const size_t n_byte);
    template <uint16_t row, uint16_t col>
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtMatrix<row, col, m_type> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtMatrix3<m_type, 3, 3> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtRotation<m_type, 3, 3> &m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const dtVector<col, m_type> &v);
    void SetRowVec(const uint16_t idxRow, const dtVector3<m_type, 3> &v);
    void SetRowVec(const uint16_t idxRow, const dtVector4<m_type, 4> &v);
    void SetRowVec(const uint16_t idxRow, const dtVector6<m_type, 6> &v);
    void SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const dtVector<row, m_type> &v);
    void SetColVec(const uint16_t idxCol, const dtVector3<m_type, 3> &v);
    void SetColVec(const uint16_t idxCol, const dtVector4<m_type, 4> &v);
    void SetColVec(const uint16_t idxCol, const dtVector6<m_type, 6> &v);
    void SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type *const GetElementsAddr() const;
    uint16_t GetRowSize() const { return m_row; } // size of row
    uint16_t GetColSize() const { return m_col; } // size of colum
    template <uint16_t row, uint16_t col>
    dtMatrix<row, col, m_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol);
    dtVector<m_col, m_type> GetRowVec(const uint16_t idxRow) const;
    dtVector<m_row, m_type> GetColVec(const uint16_t idxCol) const;
    template <uint16_t row, uint16_t col>
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, dtMatrix<row, col, m_type> &m);
    int8_t GetRowVec(const uint16_t idxRow, dtVector<m_col, m_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, dtVector<m_row, m_type> &v) const;
    dtCscMatrix<m_row, m_col, m_type> GetCscMat() const; // Compressed Sparse Column Matrix
    dtMatrix<m_col, m_row, m_type> Transpose() const;

    m_type Trace() const;
    m_type GetNorm() const;     // Frobenius Norm (Euclidean norm)
    m_type GetSqNorm() const;   // Squared Frobenius Norm (Euclidean norm)
    m_type Determinant() const; // From LU Decomposition

    dtNoPivLU<m_row, m_col, m_type> NoPivLU() const;
    dtPartialPivLU<m_row, m_col, m_type> PartialPivLU() const;
    dtLLT<m_row, m_col, m_type> LLT() const;
    dtLDLT<m_row, m_col, m_type> LDLT() const;
    dtQR<m_row, m_col, m_type> QR() const;
    dtSVD<m_row, m_col, m_type> SVD() const;

    dtMatrix<m_row, m_col, m_type> Inv(int8_t *isOk = nullptr) const;
    dtMatrix<m_col, m_row, m_type> PInv(int8_t *isOk = nullptr, m_type tolerance = std::numeric_limits<m_type>::epsilon()) const;

    /* Member access operators */
    // returns a row of modifiable elements
    inline m_type &operator()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    inline const m_type &operator()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    dtMatrix &operator=(const dtMatrix &m);                           // matrix  = matrix
    dtMatrix &operator+=(const dtMatrix &m);                          // matrix += matrix
    dtMatrix &operator-=(const dtMatrix &m);                          // matrix -= matrix
    dtMatrix &operator=(const dtMatrix3<m_type, m_row, m_col> &m);    // matrix  = matrix3
    dtMatrix &operator+=(const dtMatrix3<m_type, m_row, m_col> &m);   // matrix += matrix3
    dtMatrix &operator-=(const dtMatrix3<m_type, m_row, m_col> &m);   // matrix -= matrix3
    dtMatrix &operator=(const dtRotation<m_type, m_row, m_col> &m);   // matrix  = RotMat
    dtMatrix &operator+=(const dtRotation<m_type, m_row, m_col> &m);  // matrix += RotMat
    dtMatrix &operator-=(const dtRotation<m_type, m_row, m_col> &m);  // matrix -= RotMat
    dtMatrix &operator=(const dtTransform<m_type, m_row, m_col> &m);  // matrix + Transform
    dtMatrix &operator+=(const dtTransform<m_type, m_row, m_col> &m); // matrix + Transform
    dtMatrix &operator-=(const dtTransform<m_type, m_row, m_col> &m); // matrix - Transform
    dtMatrix &operator=(const m_type s);                              // matrix  = scalar, all elements set scalar
    dtMatrix &operator+=(const m_type s);                             // matrix += scalar, matrix(i) += scalar
    dtMatrix &operator-=(const m_type s);                             // matrix -= scalar, matrix(i) -= scalar
    dtMatrix &operator*=(const m_type s);                             // matrix *= scalar
    dtMatrix &operator/=(const m_type s);                             // matrix /= scalar
    dtCommaInit<m_row * m_col, m_type> operator<<(const m_type s);    // Init first matrix elements

    /* Arithmetic operators */
    dtMatrix operator-() const;                                           // minus sign
    dtMatrix operator+(const dtMatrix &m) const;                          // matrix + matrix
    dtMatrix operator-(const dtMatrix &m) const;                          // matrix - matrix
    dtMatrix operator+(const dtMatrix3<m_type, m_row, m_col> &m) const;   // matrix + matrix3
    dtMatrix operator-(const dtMatrix3<m_type, m_row, m_col> &m) const;   // matrix - matrix3
    dtMatrix operator+(const dtRotation<m_type, m_row, m_col> &m) const;  // matrix + RotMat
    dtMatrix operator-(const dtRotation<m_type, m_row, m_col> &m) const;  // matrix - RotMat
    dtMatrix operator+(const dtTransform<m_type, m_row, m_col> &m) const; // matrix + Transform
    dtMatrix operator-(const dtTransform<m_type, m_row, m_col> &m) const; // matrix - Transform
    dtMatrix operator+(const m_type s) const;                             // matrix + scalar, matrix(i) + scalar
    dtMatrix operator-(const m_type s) const;                             // matrix - scalar, matrix(i) - scalar
    dtMatrix operator*(const m_type s) const;                             // matrix * scalar
    dtMatrix operator/(const m_type s) const;                             // matrix / scalar

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<m_col, col, m_type> &m) const;        // matrix * matrix
    dtMatrix<m_row, m_col, m_type> operator*(const dtMatrix3<m_type, m_col, m_col> &m) const;   // matrix * matrix
    dtMatrix<m_row, m_col, m_type> operator*(const dtRotation<m_type, m_col, m_col> &m) const;  // matrix * RotMat
    dtMatrix<m_row, m_col, m_type> operator*(const dtTransform<m_type, m_col, m_col> &m) const; // matrix * Transform
    dtVector<m_row, m_type> operator*(const dtVector<m_col, m_type> &v) const;                  // matrix * vector
    dtVector<m_row, m_type> operator*(const dtVector3<m_type, m_col> &v) const;                 // matrix * vector3
    dtVector<m_row, m_type> operator*(const dtVector4<m_type, m_col> &v) const;                 // matrix * vector4
    dtVector<m_row, m_type> operator*(const dtVector6<m_type, m_col> &v) const;                 // matrix * vector6
    dtMatrix<m_row, m_col, m_type> operator&(const dtVector<m_col, m_type> &v) const;           // matrix3 * [v]x, []x is skew-symmetric matrix
    dtMatrix<m_row, m_col, m_type> operator&(const dtVector3<m_type, m_col> &v) const;          // matrix3 * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const dtMatrix &m) const;                          // (true or false) matrix1 == matrix
    bool operator!=(const dtMatrix &m) const;                          // (true or false) matrix1 != matrix
    bool operator==(const dtMatrix3<m_type, m_row, m_col> &m) const;   // (true or false) matrix == matrix3
    bool operator!=(const dtMatrix3<m_type, m_row, m_col> &m) const;   // (true or false) matrix != matrix3
    bool operator==(const dtRotation<m_type, m_row, m_col> &m) const;  // (true or false) matrix == RotMat
    bool operator!=(const dtRotation<m_type, m_row, m_col> &m) const;  // (true or false) matrix != RotMat
    bool operator==(const dtTransform<m_type, m_row, m_col> &m) const; // (true or false) matrix == Transform
    bool operator!=(const dtTransform<m_type, m_row, m_col> &m) const; // (true or false) matrix != Transform

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class dtCscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class dtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class dtTransform;

    template <uint16_t row, typename type> friend class dtVector;
    template <typename type, uint16_t row> friend class dtVector3;
    template <typename type, uint16_t row> friend class dtVector4;
    template <typename type, uint16_t row> friend class dtVector6;

    template <uint16_t row, uint16_t col, typename type> friend class dtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtQR;
    template <uint16_t row, uint16_t col, typename type> friend class dtSVD;

    /* Friend template function */
    template <uint16_t row, uint16_t col, typename type>
    friend dtMatrix<row, col, type> operator*(const type s, const dtMatrix<row, col, type> &m); // scalar * matrix
};

} // namespace dtMath

#include "dtMatrix.tpp"

#endif // DTMATH_DTMATRIX_H_