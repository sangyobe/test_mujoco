/*!
\file       dtRotation.h
\brief      dtMath, Rotation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTROTATION_H_
#define DTMATH_DTROTATION_H_

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
template <typename m_type, uint16_t m_row> class dtQuaternion;
template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class dtMatrix3;

template <typename m_type = float, uint16_t m_row = 3, uint16_t m_col = 3>
class dtRotation
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    dtRotation(const m_type *element);
    inline void Euler2RotMat(const uint16_t order, const m_type *e);
    inline void Quat2RotMat(const m_type *q);

public:
    dtRotation();
    dtRotation(const m_type *element, const size_t n_byte);
    dtRotation(const char c, const m_type *element, const size_t n_byte);
    dtRotation(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    dtRotation(const uint16_t order, const m_type angle);
    dtRotation(const uint16_t order, const m_type angle1, const m_type angle2);
    dtRotation(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    dtRotation(const dtRotation &m);
    dtRotation(const dtMatrix3<m_type, m_row, m_col> &m);
    dtRotation(const dtMatrix<m_row, m_col, m_type> &m);
    dtRotation(const uint16_t order, const dtVector3<m_type, 3> &e);
    dtRotation(const uint16_t order, const dtVector<3, m_type> &e);
    dtRotation(const dtQuaternion<m_type, 4> &q);
    ~dtRotation() {}

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

    void SetElement(const uint16_t order, const m_type angle);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);

    void SetElement(const dtRotation &m);
    void SetElement(const dtMatrix3<m_type, m_row, m_col> &m);
    void SetElement(const dtMatrix<m_row, m_col, m_type> &m);

    void SetElement(const uint16_t order, const dtVector3<m_type, 3> &e);
    void SetElement(const uint16_t order, const dtVector<3, m_type> &e);
    void SetElement(const uint16_t order, const m_type *e);

    void SetElement(const dtQuaternion<m_type, 4> &q);
    void SetElement(const m_type *q);
    void SetElement(const m_type w, const m_type x, const m_type y, const m_type z);

    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type *const GetElementsAddr() const;
    dtVector3<m_type, 3> GetRowVec(const uint16_t idxRow) const;
    dtVector3<m_type, 3> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, dtVector3<m_type, 3> &v) const;
    int8_t GetColVec(const uint16_t idxCol, dtVector3<m_type, 3> &v) const;
    dtVector3<m_type, 3> GetEulerAngles(uint16_t order) const;
    dtRotation Transpose() const;
    dtRotation log() const;                                // log(R) = ln(R) : SO(3) -> so(3), SO(3) is Special Orthogonal Group, so(3) is the set of skew-symmetric 3x3 matrices
    dtVector3<m_type, 3> Log() const;                      // Log(R) = u*phi : SO(3) -> R3
    dtRotation ode(m_type wx, m_type wy, m_type wz) const; // dR/dt = R*[w]x
    dtRotation ode(m_type *w) const;                       // dR/dt = R*[w]x
    dtRotation ode(dtVector3<m_type, 3> w) const;          // dR/dt = R*[w]x
    dtRotation ode(dtVector<3, m_type> w) const;           // dR/dt = R*[w]x
    dtRotation Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type &operator()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    const m_type &operator()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    dtRotation &operator=(const dtRotation &m);                    // matrix = matrix
    dtCommaInit<m_row * m_col, m_type> operator<<(const m_type s); // Init first matrix elements

    /* Arithmetic operators */
    dtRotation operator-() const; // minus sign
    dtMatrix3<m_type, m_row, m_col> operator+(const dtRotation &m) const;
    dtMatrix3<m_type, m_row, m_col> operator-(const dtRotation &m) const;
    dtMatrix3<m_type, m_row, m_col> operator+(const dtMatrix3<m_type, m_row, m_col> &m) const;
    dtMatrix3<m_type, m_row, m_col> operator-(const dtMatrix3<m_type, m_row, m_col> &m) const;
    dtMatrix3<m_type, m_row, m_col> operator+(const dtMatrix<m_row, m_col, m_type> &m) const;
    dtMatrix3<m_type, m_row, m_col> operator-(const dtMatrix<m_row, m_col, m_type> &m) const;
    dtMatrix3<m_type, m_row, m_col> operator*(const m_type s) const;
    dtMatrix3<m_type, m_row, m_col> operator/(const m_type s) const;

    template <uint16_t col>
    dtMatrix<m_row, col, m_type> operator*(const dtMatrix<m_col, col, m_type> &m) const;       // RotMat * matrix
    dtMatrix3<m_type, m_row, m_col> operator*(const dtMatrix3<m_type, m_row, m_col> &m) const; // RotMat * matrix
    dtRotation operator*(const dtRotation &m) const;                                           // RotMat * RotMat
    dtVector<m_row, m_type> operator*(const dtVector<m_col, m_type> &v) const;                 // RotMat * vector
    dtVector3<m_type, m_row> operator*(const dtVector3<m_type, m_col> &v) const;               // RotMat * vector
    dtRotation operator&(const dtVector<m_col, m_type> &v) const;                              // RotMat * [v]x, []x is skew-symmetric matrix
    dtRotation operator&(const dtVector3<m_type, m_col> &v) const;                             // RotMat * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const dtRotation &m) const;                      // (true or false) matrix == matrix
    bool operator!=(const dtRotation &m) const;                      // (true or false) matrix != matrix
    bool operator==(const dtMatrix3<m_type, m_row, m_col> &m) const; // (true or false) matrix == matrix
    bool operator!=(const dtMatrix3<m_type, m_row, m_col> &m) const; // (true or false) matrix != matrix
    bool operator==(const dtMatrix<m_row, m_col, m_type> &m) const;  // (true or false) matrix == matrix
    bool operator!=(const dtMatrix<m_row, m_col, m_type> &m) const;  // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class dtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class dtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class dtTransform;
    template <uint16_t row, uint16_t col, typename type> friend class dtMatrix;

    template <typename type, uint16_t row> friend class dtQuaternion;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend dtMatrix3<type, row, col> operator*(const type s, const dtRotation<type, row, col> &m); // scalar * RotMat
};

} // namespace dtMath

#include "dtRotation.tpp"

#endif // DTMATH_DTROTATION_H_
