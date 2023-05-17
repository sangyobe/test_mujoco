/*!
\file       dtVector.h
\brief      dtMath, General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR_H_
#define DTMATH_DTVECTOR_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
// #include <cstdarg>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_size, typename m_type> class dtCommaInit;
template <typename m_type, uint16_t m_row> class dtVector3;
template <typename m_type, uint16_t m_row> class dtVector4;
template <typename m_type, uint16_t m_row> class dtVector6;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtCscMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtRotation;

template <uint16_t m_row, typename m_type = float>
class dtVector
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    dtVector(const m_type *element);
    inline void CrossProduct(const m_type *v);

public:
    dtVector();
    dtVector(const m_type *element, const size_t n_byte);
    dtVector(const dtVector &v);
    dtVector(const dtMatrix<m_row, 1, m_type> &v);
    ~dtVector() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type *element, const size_t n_byte);
    void SetElement(const dtVector &v);
    void SetElement(const dtMatrix<m_row, 1, m_type> &v);
    // void SetElement(const m_type elem0, ...);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const dtVector<row, m_type> &v);
    void SetBlock(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    void SetBlock(const uint16_t idxRow, const dtVector3<m_type, 3> &v);
    void SetBlock(const uint16_t idxRow, const dtVector4<m_type, 4> &v);
    void SetBlock(const uint16_t idxRow, const dtVector6<m_type, 6> &v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const dtMatrix<row, 1, m_type> &v);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const m_type *const GetElementsAddr() const;
    uint16_t GetDim() const { return m_row; }
    template <uint16_t row>
    dtVector<row, m_type> GetBlock(const uint16_t idx);
    dtVector3<m_type, 3> GetBlockVec3(const uint16_t idx);
    dtVector4<m_type, 4> GetBlockVec4(const uint16_t idx);
    dtVector6<m_type, 6> GetBlockVec6(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, dtVector<row, m_type> &v);
    int8_t GetBlockVec3(const uint16_t idx, dtVector3<m_type, 3> &v);
    int8_t GetBlockVec4(const uint16_t idx, dtVector4<m_type, 4> &v);
    int8_t GetBlockVec6(const uint16_t idx, dtVector6<m_type, 6> &v);
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    dtVector GetNormalized() const;
    dtMatrix<3, 3, m_type> GetSkew() const;
    dtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    dtVector &operator=(const dtVector &v);                    // vector1  = vector2
    dtVector &operator+=(const dtVector &v);                   // vector1 += vector2
    dtVector &operator-=(const dtVector &v);                   // vector1 -= vector2
    dtVector &operator*=(const dtVector &v);                   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    dtVector &operator/=(const dtVector &v);                   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    dtVector &operator=(const dtVector3<m_type, m_row> &v);    // vector1  = vector2
    dtVector &operator+=(const dtVector3<m_type, m_row> &v);   // vector1 += vector2
    dtVector &operator-=(const dtVector3<m_type, m_row> &v);   // vector1 -= vector2
    dtVector &operator*=(const dtVector3<m_type, m_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    dtVector &operator/=(const dtVector3<m_type, m_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    dtVector &operator=(const dtVector4<m_type, m_row> &v);    // vector1  = vector2
    dtVector &operator+=(const dtVector4<m_type, m_row> &v);   // vector1 += vector2
    dtVector &operator-=(const dtVector4<m_type, m_row> &v);   // vector1 -= vector2
    dtVector &operator*=(const dtVector4<m_type, m_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    dtVector &operator/=(const dtVector4<m_type, m_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    dtVector &operator=(const dtVector6<m_type, m_row> &v);    // vector1  = vector2
    dtVector &operator+=(const dtVector6<m_type, m_row> &v);   // vector1 += vector2
    dtVector &operator-=(const dtVector6<m_type, m_row> &v);   // vector1 -= vector2
    dtVector &operator*=(const dtVector6<m_type, m_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    dtVector &operator/=(const dtVector6<m_type, m_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    dtVector &operator=(const dtMatrix<m_row, 1, m_type> &v);  // vector1  = vector2
    dtVector &operator+=(const dtMatrix<m_row, 1, m_type> &v); // vector1 += vector2
    dtVector &operator-=(const dtMatrix<m_row, 1, m_type> &v); // vector1 -= vector2
    dtVector &operator*=(const dtMatrix<m_row, 1, m_type> &v); // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    dtVector &operator/=(const dtMatrix<m_row, 1, m_type> &v); // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    dtVector &operator=(const m_type s);                       // vector1  = scalar, all elements set scalar
    dtVector &operator+=(const m_type s);                      // vector1 += scalar, vector1(i) = vector1(i) + scalar
    dtVector &operator-=(const m_type s);                      // vector1 -= scalar, vector1(i) = vector1(i) - scalar
    dtVector &operator*=(const m_type s);                      // vector1 *= scalar
    dtVector &operator/=(const m_type s);                      // vector1 /= scalar
    dtVector &operator&=(const dtVector3<m_type, m_row> &v);   // vector1 x vector2 cross product
    dtVector &operator&=(const dtVector<m_row, m_type> &v);    // vector1 x vector2 cross product
    dtVector &operator&=(const dtMatrix<m_row, 1, m_type> &v); // vector1 x vector2 cross product
    dtCommaInit<m_row, m_type> operator<<(const m_type s);     // Init first matrix elements

    /* Arithmetic operators */
    dtVector operator-() const;                                                                  // minus sign
    dtVector operator+(const dtVector &v) const;                                                 // vector + vector
    dtVector operator-(const dtVector &v) const;                                                 // vector - vector
    dtVector operator*(const dtVector &v) const;                                                 // vector * vector = m_elem(i) * v.m_elem(i)
    dtVector operator/(const dtVector &v) const;                                                 // vector / vector = m_elem(i) / v.m_elem(i)
    dtVector operator+(const dtVector3<m_type, m_row> &v) const;                                 // vector + vector
    dtVector operator-(const dtVector3<m_type, m_row> &v) const;                                 // vector - vector
    dtVector operator*(const dtVector3<m_type, m_row> &v) const;                                 // vector * vector = m_elem(i) * v.m_elem(i)
    dtVector operator/(const dtVector3<m_type, m_row> &v) const;                                 // vector / vector = m_elem(i) / v.m_elem(i)
    dtVector operator+(const dtVector4<m_type, m_row> &v) const;                                 // vector + vector
    dtVector operator-(const dtVector4<m_type, m_row> &v) const;                                 // vector - vector
    dtVector operator*(const dtVector4<m_type, m_row> &v) const;                                 // vector * vector = m_elem(i) * v.m_elem(i)
    dtVector operator/(const dtVector4<m_type, m_row> &v) const;                                 // vector / vector = m_elem(i) / v.m_elem(i)
    dtVector operator+(const dtVector6<m_type, m_row> &v) const;                                 // vector + vector
    dtVector operator-(const dtVector6<m_type, m_row> &v) const;                                 // vector - vector
    dtVector operator*(const dtVector6<m_type, m_row> &v) const;                                 // vector * vector = m_elem(i) * v.m_elem(i)
    dtVector operator/(const dtVector6<m_type, m_row> &v) const;                                 // vector / vector = m_elem(i) / v.m_elem(i)
    dtVector operator+(const dtMatrix<m_row, 1, m_type> &v) const;                               // vector + vector
    dtVector operator-(const dtMatrix<m_row, 1, m_type> &v) const;                               // vector - vector
    dtVector operator*(const dtMatrix<m_row, 1, m_type> &v) const;                               // vector * vector = m_elem(i) * v.m_elem(i)
    dtVector operator/(const dtMatrix<m_row, 1, m_type> &v) const;                               // vector / vector = m_elem(i) / v.m_elem(i)
    dtVector operator+(const m_type s) const;                                                    // vector + scalar = m_elem(i) + scalar
    dtVector operator-(const m_type s) const;                                                    // vector - scalar = m_elem(i) - scalar
    dtVector operator*(const m_type s) const;                                                    // vector * scalar
    dtVector operator/(const m_type s) const;                                                    // vector / scalar
    dtVector operator&(const dtVector3<m_type, m_row> &v) const;                                 // vector x vector cross product
    dtVector operator&(const dtVector<m_row, m_type> &v) const;                                  // vector x vector cross product
    dtVector operator&(const dtMatrix<m_row, 1, m_type> &v) const;                               // vector x vector cross product
    dtMatrix3<m_type, m_row, m_row> operator&(const dtMatrix3<m_type, m_row, m_row> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    dtRotation<m_type, m_row, m_row> operator&(const dtRotation<m_type, m_row, m_row> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<1, col, m_type> &m) const; // vector1 * matrix(1xcol) outer product

    m_type dot(const dtVector &v) const;                   // vector1 * vector2 dot(inner) product
    m_type dot(const dtVector3<m_type, m_row> &v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const dtVector4<m_type, m_row> &v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const dtVector6<m_type, m_row> &v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const dtMatrix<m_row, 1, m_type> &v) const; // vector1 * vector2 dot(inner) product

    /* Comparison operators */
    bool operator==(const dtVector &v) const;                   // (true or false) vector1 == vector2
    bool operator!=(const dtVector &v) const;                   // (true or false) vector1 != vector2
    bool operator==(const dtVector3<m_type, m_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const dtVector3<m_type, m_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const dtVector4<m_type, m_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const dtVector4<m_type, m_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const dtVector6<m_type, m_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const dtVector6<m_type, m_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, typename type> friend class dtVector;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class dtCscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class dtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;
    template <typename type, uint16_t row> friend class dtVector3;
    template <typename type, uint16_t row> friend class dtVector4;
    template <typename type, uint16_t row> friend class dtVector6;
    template <typename type, uint16_t row> friend class dtQuaternion;

    template <uint16_t row, uint16_t col, typename type> friend class dtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class dtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class dtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class dtSVD;

    /* Friend template function */
    template <uint16_t row, typename type>
    friend dtVector<row, type> operator+(const type s, const dtVector<row, type> &v); // scalar + vector = scalar + m_elem(i)
    template <uint16_t row, typename type>
    friend dtVector<row, type> operator-(const type s, const dtVector<row, type> &v); // scalar - vector = scalar - m_elem(i)
    template <uint16_t row, typename type>
    friend dtVector<row, type> operator*(const type s, const dtVector<row, type> &v); // scalar * vector = scalar * m_elem(i)
    template <uint16_t row, typename type>
    friend dtVector<row, type> operator/(const type s, const dtVector<row, type> &v); // scalar / vector = scalar / m_elem(i)
};

} // namespace dtMath

#include "dtVector.tpp"

#endif // DTMATH_DTVECTOR_H_
