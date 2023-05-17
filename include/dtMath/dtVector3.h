/*!
\file       dtVector3.h
\brief      dtMath, 3x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR3_H_
#define DTMATH_DTVECTOR3_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_size, typename m_type> class dtCommaInit;
template <uint16_t m_row, typename m_type> class dtVector;
template <typename m_type, uint16_t m_row> class dtVector4;
template <typename m_type, uint16_t m_row> class dtVector6;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtRotation;

template <typename m_type = float, uint16_t m_row = 3>
class dtVector3
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    dtVector3(const m_type *element);
    inline void CrossProduct(const m_type *v);

public:
    dtVector3();
    dtVector3(const m_type *element, const size_t n_byte);
    dtVector3(const m_type i, const m_type j, const m_type k);
    dtVector3(const dtVector3 &v);
    dtVector3(const dtVector<m_row, m_type> &v);
    dtVector3(const dtMatrix<m_row, 1, m_type> &v);
    ~dtVector3() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type *element, const size_t n_byte);
    void SetElement(const m_type i, const m_type j, const m_type k);
    void SetElement(const dtVector3 &v);
    void SetElement(const dtVector<m_row, m_type> &v);
    void SetElement(const dtMatrix<m_row, 1, m_type> &v);
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
    template <uint16_t row>
    dtVector<row, m_type> GetBlock(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, dtVector<row, m_type> &v);
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    dtVector3 GetNormalized() const;
    dtMatrix3<m_type, m_row, m_row> GetSkew() const;
    dtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    dtVector3 &operator=(const dtVector3 &v);                   // vector1  = vector2
    dtVector3 &operator+=(const dtVector3 &v);                  // vector1 += vector2
    dtVector3 &operator-=(const dtVector3 &v);                  // vector1 -= vector2
    dtVector3 &operator*=(const dtVector3 &v);                  // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector3 &operator/=(const dtVector3 &v);                  // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector3 &operator=(const dtVector<m_row, m_type> &v);     // vector1  = vector2
    dtVector3 &operator+=(const dtVector<m_row, m_type> &v);    // vector1 += vector2
    dtVector3 &operator-=(const dtVector<m_row, m_type> &v);    // vector1 -= vector2
    dtVector3 &operator*=(const dtVector<m_row, m_type> &v);    // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector3 &operator/=(const dtVector<m_row, m_type> &v);    // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector3 &operator=(const dtMatrix<m_row, 1, m_type> &v);  // vector1  = vector2
    dtVector3 &operator+=(const dtMatrix<m_row, 1, m_type> &v); // vector1 += vector2
    dtVector3 &operator-=(const dtMatrix<m_row, 1, m_type> &v); // vector1 -= vector2
    dtVector3 &operator*=(const dtMatrix<m_row, 1, m_type> &v); // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector3 &operator/=(const dtMatrix<m_row, 1, m_type> &v); // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector3 &operator=(const m_type s);                       // vector1  = scalar, all elements set scalar
    dtVector3 &operator+=(const m_type s);                      // vector1 += scalar, vector1(i) += scalar
    dtVector3 &operator-=(const m_type s);                      // vector1 -= scalar, vector1(i) -= scalar
    dtVector3 &operator*=(const m_type s);                      // vector1 *= scalar
    dtVector3 &operator/=(const m_type s);                      // vector1 /= scalar
    dtVector3 &operator&=(const dtVector3 &v);                  // vector1 x vector2 cross product
    dtVector3 &operator&=(const dtVector<m_row, m_type> &v);    // vector1 x vector2 cross product
    dtVector3 &operator&=(const dtMatrix<m_row, 1, m_type> &v); // vector1 x vector2 cross product
    dtCommaInit<m_row, m_type> operator<<(const m_type s);      // Init first matrix elements

    /* Arithmetic operators */
    dtVector3 operator-() const;                                                                 // minus sign
    dtVector3 operator+(const dtVector3 &v) const;                                               // vector + vector
    dtVector3 operator-(const dtVector3 &v) const;                                               // vector - vector
    dtVector3 operator*(const dtVector3 &v) const;                                               // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector3 operator/(const dtVector3 &v) const;                                               // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector3 operator+(const dtVector<m_row, m_type> &v) const;                                 // vector + vector
    dtVector3 operator-(const dtVector<m_row, m_type> &v) const;                                 // vector - vector
    dtVector3 operator*(const dtVector<m_row, m_type> &v) const;                                 // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector3 operator/(const dtVector<m_row, m_type> &v) const;                                 // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector3 operator+(const dtMatrix<m_row, 1, m_type> &v) const;                              // vector + vector
    dtVector3 operator-(const dtMatrix<m_row, 1, m_type> &v) const;                              // vector - vector
    dtVector3 operator*(const dtMatrix<m_row, 1, m_type> &v) const;                              // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector3 operator/(const dtMatrix<m_row, 1, m_type> &v) const;                              // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector3 operator+(const m_type s) const;                                                   // vector + scalar, vector(i) = vector1(i) + scalar
    dtVector3 operator-(const m_type s) const;                                                   // vector - scalar, vector(i) = vector1(i) - scalar
    dtVector3 operator*(const m_type s) const;                                                   // vector * scalar
    dtVector3 operator/(const m_type s) const;                                                   // vector / scalar
    dtVector3 operator&(const dtVector3 &v) const;                                               // vector x vector cross product
    dtVector3 operator&(const dtVector<m_row, m_type> &v) const;                                 // vector x vector cross product
    dtVector3 operator&(const dtMatrix<m_row, 1, m_type> &v) const;                              // vector x vector cross product
    dtMatrix3<m_type, m_row, m_row> operator&(const dtMatrix3<m_type, m_row, m_row> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    dtRotation<m_type, m_row, m_row> operator&(const dtRotation<m_type, m_row, m_row> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<1, col, m_type> &m) const; // vector1 * matrix(1xcol) outer product

    m_type dot(const dtVector3 &v) const;                  // vector1 * vector2 dot(inner) product
    m_type dot(const dtVector<m_row, m_type> &v) const;    // vector1 * vector2 dot(inner) product
    m_type dot(const dtMatrix<m_row, 1, m_type> &v) const; // vector1 * vector2 dot(inner) product

    /* Comparison operators */
    bool operator==(const dtVector3 &v) const;                  // (true or false) vector1 == vector2
    bool operator!=(const dtVector3 &v) const;                  // (true or false) vector1 != vector2
    bool operator==(const dtVector<m_row, m_type> &v) const;    // (true or false) vector1 == vector2
    bool operator!=(const dtVector<m_row, m_type> &v) const;    // (true or false) vector1 != vector2
    bool operator==(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class dtVector3;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class dtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class dtTransform;

    template <uint16_t row, typename type> friend class dtVector;
    template <typename type, uint16_t row> friend class dtVector4;
    template <typename type, uint16_t row> friend class dtVector6;
    template <typename type, uint16_t row> friend class dtQuaternion;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend dtVector3<type, row> operator+(const type s, const dtVector3<type, row> &v); // scalar + vector, scalar + vector(i)
    template <typename type, uint16_t row>
    friend dtVector3<type, row> operator-(const type s, const dtVector3<type, row> &v); // scalar - vector, scalar - vector(i)
    template <typename type, uint16_t row>
    friend dtVector3<type, row> operator*(const type s, const dtVector3<type, row> &v); // scalar * vector, scalar * vector(i)
    template <typename type, uint16_t row>
    friend dtVector3<type, row> operator/(const type s, const dtVector3<type, row> &v); // scalar / vector, scalar / vector(i)
};

} // namespace dtMath

#include "dtVector3.tpp"

#endif // DTMATH_DTVECTOR3_H_
