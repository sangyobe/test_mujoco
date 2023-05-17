/*!
\file       dtTransform.h
\brief      dtMath, Homogeneous transformation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTTRANSFORM_H_
#define DTMATH_DTTRANSFORM_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>

namespace dtMath
{

template <uint16_t m_size, typename m_type> class dtCommaInit;
template <uint16_t m_row, typename m_type> class dtVector;
template <typename m_type, uint16_t m_row> class dtVector3;
template <typename m_type, uint16_t m_row> class dtVector6;
template <typename m_type, uint16_t m_row> class dtQuaternion;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;

template <typename m_type, uint16_t m_row, uint16_t m_col>
class dtRotation;

template <typename m_type = float, uint16_t m_row = 4, uint16_t m_col = 4>
class dtTransform
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    dtRotation<m_type, 3, 3> m_R;
    dtVector3<m_type, 3> m_p;
    m_type m_dummy[4];

public:
    dtTransform();
    dtTransform(const dtRotation<m_type, 3, 3> &R, const dtVector3<m_type, 3> &p);
    dtTransform(const dtQuaternion<m_type, 4> &q, const dtVector3<m_type, 3> &p);
    dtTransform(const uint16_t order, const dtVector3<m_type, 3> &e, const dtVector3<m_type, 3> &p);
    dtTransform(const dtTransform &m);
    ~dtTransform() {}

    void SetZero();
    void SetIdentity();
    void SetElement(const dtVector3<m_type, 3> &p);
    void SetElement(const dtRotation<m_type, 3, 3> &R);
    void SetElement(const dtQuaternion<m_type, 4> &q);
    void SetElement(const uint16_t order, const dtVector3<m_type, 3> &e);
    void SetElement(const dtRotation<m_type, 3, 3> &R, const dtVector3<m_type, 3> &p);
    void SetElement(const dtQuaternion<m_type, 4> &q, const dtVector3<m_type, 3> &p);
    void SetElement(const uint16_t order, const dtVector3<m_type, 3> &e, const dtVector3<m_type, 3> &p);
    void SetElement(const dtTransform &m);

    dtQuaternion<m_type, 4> q() const { return dtQuaternion<m_type, 4>(m_R); }
    dtRotation<m_type, 3, 3> R() const { return m_R; }
    dtVector3<m_type, 3> e(uint16_t order) const { return m_R.GetEulerAngles(order); }
    dtVector3<m_type, 3> p() const { return m_p; }
    dtVector6<m_type, 6> GetError(const dtTransform &m) const;
    dtMatrix<m_col, m_row, m_type> Transpose() const;
    dtTransform Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow, uint16_t icol);
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow, uint16_t icol) const;

    /* Assignment operators */
    dtTransform &operator=(const dtTransform &m); // transform matrix  = transform matrix

    /* Arithmetic operators */
    dtMatrix<m_row, m_col, m_type> operator+(const dtMatrix<m_row, m_col, m_type> &m) const; // transform matrix + matrix
    dtMatrix<m_row, m_col, m_type> operator-(const dtMatrix<m_row, m_col, m_type> &m) const; // transform matrix - matrix

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<m_row, col, m_type> &m) const; // transform matrix * matrix
    dtTransform operator*(const dtTransform &m) const;                                   // transform matrix * transform matrix
    dtVector3<m_type, 3> operator*(const dtVector<3, m_type> &v) const;                  // matrix * vector
    dtVector3<m_type, 3> operator*(const dtVector3<m_type, 3> &v) const;                 // matrix * vector3

    /* Comparison operators */
    bool operator==(const dtTransform &m) const; // (true or false) matrix1 == matrix
    bool operator!=(const dtTransform &m) const; // (true or false) matrix1 == matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class dtTransform;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;
};

} // namespace dtMath

#include "dtTransform.tpp"

#endif // DTMATH_DTTRANSFORM_H_
