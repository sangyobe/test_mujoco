/*!
\file       dtVector6.h
\brief      dtMath, 6x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR6_H_
#define DTMATH_DTVECTOR6_H_

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
template <typename m_type, uint16_t m_row> class dtVector3;
template <typename m_type, uint16_t m_row> class dtVector4;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;

template <typename m_type = float, uint16_t m_row = 6>
class dtVector6
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    dtVector6(const m_type *element);

public:
    dtVector6();
    dtVector6(const m_type *element, const size_t n_byte);
    dtVector6(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz);
    dtVector6(const dtVector6 &v);
    dtVector6(const dtVector3<m_type, 3> &p, const dtVector3<m_type, 3> &o);
    dtVector6(const dtVector<m_row, m_type> &v);
    dtVector6(const dtMatrix<m_row, 1, m_type> &v);
    ~dtVector6() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type *element, const size_t n_byte);
    void SetElement(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz);
    void SetElement(const dtVector6 &v);
    void SetElement(const dtVector3<m_type, 3> &p, const dtVector3<m_type, 3> &o);
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
    dtVector3<m_type, 3> GetBlockVec3(const uint16_t idx);
    dtVector4<m_type, 4> GetBlockVec4(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, dtVector<row, m_type> &v);
    int8_t GetBlockVec3(const uint16_t idx, dtVector3<m_type, 3> &v);
    int8_t GetBlockVec4(const uint16_t idx, dtVector4<m_type, 4> &v);
    dtVector3<m_type, 3> GetPos() const;
    dtVector3<m_type, 3> GetOri() const;
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    dtVector6 GetNormalized() const;
    dtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    dtVector6 &operator=(const dtVector6 &v);                   // vector1  = vector2
    dtVector6 &operator+=(const dtVector6 &v);                  // vector1 += vector2
    dtVector6 &operator-=(const dtVector6 &v);                  // vector1 -= vector2
    dtVector6 &operator*=(const dtVector6 &v);                  // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector6 &operator/=(const dtVector6 &v);                  // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector6 &operator=(const dtVector<m_row, m_type> &v);     // vector1  = vector2
    dtVector6 &operator+=(const dtVector<m_row, m_type> &v);    // vector1 += vector2
    dtVector6 &operator-=(const dtVector<m_row, m_type> &v);    // vector1 -= vector2
    dtVector6 &operator*=(const dtVector<m_row, m_type> &v);    // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector6 &operator/=(const dtVector<m_row, m_type> &v);    // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector6 &operator=(const dtMatrix<m_row, 1, m_type> &v);  // vector1  = vector2
    dtVector6 &operator+=(const dtMatrix<m_row, 1, m_type> &v); // vector1 += vector2
    dtVector6 &operator-=(const dtMatrix<m_row, 1, m_type> &v); // vector1 -= vector2
    dtVector6 &operator*=(const dtMatrix<m_row, 1, m_type> &v); // vector1 *= vector2, vector1(i) *= vector2(i)
    dtVector6 &operator/=(const dtMatrix<m_row, 1, m_type> &v); // vector1 /= vector2, vector1(i) /= vector2(i)
    dtVector6 &operator=(const m_type s);                       // vector1  = scalar, all elements set scalar
    dtVector6 &operator+=(const m_type s);                      // vector1 += scalar, vector1(i) += scalar
    dtVector6 &operator-=(const m_type s);                      // vector1 -= scalar, vector1(i) -= scalar
    dtVector6 &operator*=(const m_type s);                      // vector1 *= scalar
    dtVector6 &operator/=(const m_type s);                      // vector1 /= scalar
    dtCommaInit<m_row, m_type> operator<<(const m_type s);      // Init first matrix elements

    /* Arithmetic operators */
    dtVector6 operator-() const;                                    // minus sign
    dtVector6 operator+(const dtVector6 &v) const;                  // vector + vector
    dtVector6 operator-(const dtVector6 &v) const;                  // vector - vector
    dtVector6 operator*(const dtVector6 &v) const;                  // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector6 operator/(const dtVector6 &v) const;                  // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector6 operator+(const dtVector<m_row, m_type> &v) const;    // vector + vector
    dtVector6 operator-(const dtVector<m_row, m_type> &v) const;    // vector - vector
    dtVector6 operator*(const dtVector<m_row, m_type> &v) const;    // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector6 operator/(const dtVector<m_row, m_type> &v) const;    // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector6 operator+(const dtMatrix<m_row, 1, m_type> &v) const; // vector + vector
    dtVector6 operator-(const dtMatrix<m_row, 1, m_type> &v) const; // vector - vector
    dtVector6 operator*(const dtMatrix<m_row, 1, m_type> &v) const; // vector * vector, vector(i) = vector1(i) * vector2(i)
    dtVector6 operator/(const dtMatrix<m_row, 1, m_type> &v) const; // vector / vector, vector(i) = vector1(i) / vector2(i)
    dtVector6 operator+(const m_type s) const;                      // vector + scalar, vector(i) = vector1(i) + scalar
    dtVector6 operator-(const m_type s) const;                      // vector - scalar, vector(i) = vector1(i) - scalar
    dtVector6 operator*(const m_type s) const;                      // vector * scalar
    dtVector6 operator/(const m_type s) const;                      // vector / scalar

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<1, col, m_type> &m) const;

    m_type dot(const dtVector6 &v) const;                  // vector * vector dot(inner) product
    m_type dot(const dtVector<m_row, m_type> &v) const;    // vector * vector dot(inner) product
    m_type dot(const dtMatrix<m_row, 1, m_type> &v) const; // Vector * vector dot(inner) product

    /* Comparison operators */
    bool operator==(const dtVector6 &v) const;                  // (true or false) vector1 == vector2
    bool operator!=(const dtVector6 &v) const;                  // (true or false) vector1 != vector2
    bool operator==(const dtVector<m_row, m_type> &v) const;    // (true or false) vector1 == vector2
    bool operator!=(const dtVector<m_row, m_type> &v) const;    // (true or false) vector1 != vector2
    bool operator==(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const dtMatrix<m_row, 1, m_type> &v) const; // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class dtVector6;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
    template <uint16_t row, typename type> friend class dtVector;
    template <typename type, uint16_t row> friend class dtVector3;
    template <typename type, uint16_t row> friend class dtVector4;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend dtVector6<type, row> operator+(const type s, const dtVector6<type, row> &v); // scalar + vector, scalar + vector(i)
    template <typename type, uint16_t row>
    friend dtVector6<type, row> operator-(const type s, const dtVector6<type, row> &v); // scalar - vector, scalar - vector(i)
    template <typename type, uint16_t row>
    friend dtVector6<type, row> operator*(const type s, const dtVector6<type, row> &v); // scalar * vector, scalar * vector(i)
    template <typename type, uint16_t row>
    friend dtVector6<type, row> operator/(const type s, const dtVector6<type, row> &v); // scalar / vector, scalar / vector(i)
};

} // namespace dtMath

#include "dtVector6.tpp"

#endif // DTMATH_DTVECTOR6_H_
