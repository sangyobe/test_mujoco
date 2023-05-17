/*!
\file       dtQuaternion.h
\brief      dtMath, Quaternion class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQUATERNION_H_
#define DTMATH_DTQUATERNION_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
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
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtRotation;

template <typename m_type = float, uint16_t m_row = 4>
class dtQuaternion
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row]; // [w, x, y, z]
    inline int SGN(const m_type x) const { return (x >= m_type(0)) ? 1 : -1; }
    inline void Euler2Quat(const uint16_t order, const m_type *e);
    inline void RotMat2Quat(const m_type *rm);

public:
    dtQuaternion();
    dtQuaternion(const m_type *element);
    dtQuaternion(const m_type w, const m_type x, const m_type y, const m_type z);
    dtQuaternion(const uint16_t order, const m_type angle);
    dtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2);
    dtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    dtQuaternion(const dtQuaternion &q);
    dtQuaternion(const uint16_t order, const dtVector3<m_type, 3> &e);
    dtQuaternion(const uint16_t order, const dtVector<3, m_type> &e);
    dtQuaternion(const dtRotation<m_type, 3, 3> &rm);
    ~dtQuaternion() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type *element);
    void SetElement(const m_type w, const m_type x, const m_type y, const m_type z);
    void SetElement(const uint16_t order, const m_type angle);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    void SetElement(const dtQuaternion &q);
    void SetElement(const uint16_t order, const dtVector3<m_type, 3> &e);
    void SetElement(const uint16_t order, const dtVector<3, m_type> &e);
    void SetElement(const dtRotation<m_type, 3, 3> &rm);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const m_type *const GetElementsAddr() const;
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    dtQuaternion GetNormalized() const;
    dtQuaternion GetConj() const;
    dtVector3<m_type, 3> GetEulerAngles(const uint16_t order) const;
    dtVector3<m_type, 3> GetOriErr(const dtQuaternion &q) const;
    dtQuaternion exp() const;                                // exp(q) = e^(q), Expoential of general quaternions
    dtQuaternion log() const;                                // log(q) = ln(q), Logarithm of general quaternions
    dtVector3<m_type, 3> Log() const;                        // Log(q) = u*phi : S3 -> R3, here S3 is the 3-dimensional surface of the unit sphere of R4 = unit quaternions
    dtQuaternion ode(m_type wx, m_type wy, m_type wz) const; // dq/dt
    dtQuaternion ode(m_type *w) const;                       // dq/dt
    dtQuaternion ode(dtVector3<m_type, 3> w) const;          // dq/dt
    dtQuaternion ode(dtVector<3, m_type> w) const;           // dq/dt
    dtQuaternion Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    dtQuaternion &operator=(const dtQuaternion &q);        // quaternion1  = quaternion2
    dtQuaternion &operator+=(const dtQuaternion &q);       // quaternion1 += quaternion2
    dtQuaternion &operator-=(const dtQuaternion &q);       // quaternion1 -= quaternion2
    dtQuaternion &operator*=(const m_type s);              // quaternion1 *= scalar
    dtQuaternion &operator/=(const m_type s);              // quaternion1 /= scalar
    dtCommaInit<m_row, m_type> operator<<(const m_type s); // Init first matrix elements

    /* Arithmetic operators */
    dtQuaternion operator-() const;                      // -quaternion : minus sign
    dtQuaternion operator+(const dtQuaternion &q) const; // quaternion + quaternion
    dtQuaternion operator-(const dtQuaternion &q) const; // quaternion - quaternion
    dtQuaternion operator*(const m_type s) const;        // quaternion * scalar
    dtQuaternion operator/(const m_type s) const;        // quaternion / scalar
    dtQuaternion operator*(const dtQuaternion &q) const; // quaternion1 * quaternion2 : Quaternion Multiplication

    /* Comparison operators */
    bool operator==(const dtQuaternion &q) const; // (true or false) quaternion1 == quaternion2
    bool operator!=(const dtQuaternion &q) const; // (true or false) quaternion1 != quaternion2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class dtQuaternion;
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend dtQuaternion<type, row> operator*(const type s, const dtQuaternion<type, row> &v); // scalar * quaternion
};

} // namespace dtMath

#include "dtQuaternion.tpp"

#endif // DTMATH_DTQUATERNION_H_
