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

#ifndef DTMATH_DTROTATION_TPP_
#define DTMATH_DTROTATION_TPP_

#include "dtRotation.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation()
{
    m_elem[0] = 1;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 1;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const m_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
    m_elem[4] = element[4];
    m_elem[5] = element[5];
    m_elem[6] = element[6];
    m_elem[7] = element[7];
    m_elem[8] = element[8];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const m_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(m_type) * m_row * m_col;

    if (matSz <= n_byte)
    {
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        m_elem[6] = element[6];
        m_elem[7] = element[7];
        m_elem[8] = element[8];
    }
    else
    {
        memset(m_elem, 0, matSz);
        memcpy(m_elem, element, n_byte);
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(char c, const m_type *element, const size_t n_byte)
{
    if (c == 'a')
    {
        size_t matSz = sizeof(m_type) * m_row * m_col;

        if (matSz <= n_byte)
        {
            m_elem[0] = element[0];
            m_elem[1] = element[1];
            m_elem[2] = element[2];
            m_elem[3] = element[3];
            m_elem[4] = element[4];
            m_elem[5] = element[5];
            m_elem[6] = element[6];
            m_elem[7] = element[7];
            m_elem[8] = element[8];
        }
        else
        {
            memset(m_elem, 0, matSz);
            memcpy(m_elem, element, n_byte);
        }
    }

    else if (c == 'd')
    {
        switch (n_byte / sizeof(m_type))
        {
        case 1:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = 0;
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = 0;
            break;
        case 2:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = element[1];
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = 0;
            break;
        default:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = element[1];
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = element[2];
            break;
        }
    }

    else
    {
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(
    const m_type m00, const m_type m01, const m_type m02,
    const m_type m10, const m_type m11, const m_type m12,
    const m_type m20, const m_type m21, const m_type m22)
{
    m_elem[0] = m00;
    m_elem[1] = m01;
    m_elem[2] = m02;
    m_elem[3] = m10;
    m_elem[4] = m11;
    m_elem[5] = m12;
    m_elem[6] = m20;
    m_elem[7] = m21;
    m_elem[8] = m22;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const uint16_t order, const m_type angle)
{
    m_type c = std::cos(angle);
    m_type s = std::sin(angle);

    switch (order)
    {
    case 0x0: // X-axis
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = c;
        m_elem[5] = -s;
        m_elem[6] = 0;
        m_elem[7] = s;
        m_elem[8] = c;
        break;
    case 0x1: // Y-axis
        m_elem[0] = c;
        m_elem[1] = 0;
        m_elem[2] = s;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = -s;
        m_elem[7] = 0;
        m_elem[8] = c;
        break;
    case 0x2: // Z-axis
        m_elem[0] = c;
        m_elem[1] = -s;
        m_elem[2] = 0;
        m_elem[3] = s;
        m_elem[4] = c;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
        break;
    default:
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const uint16_t order, const m_type angle1, const m_type angle2)
{
    SetElement(order & 0xF, angle1); // R1
    dtRotation<m_type> R2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * R2;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3)
{
    SetElement(order & 0xF, angle1); // R1
    dtRotation<m_type> R2((order >> 4) & 0xF, angle2);
    dtRotation<m_type> R3((order >> 8) & 0xF, angle3);
    (*this) = (*this) * R2 * R3;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const dtRotation &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const dtMatrix3<m_type, m_row, m_col> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const dtMatrix<m_row, m_col, m_type> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const uint16_t order, const dtVector3<m_type, 3> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const uint16_t order, const dtVector<3, m_type> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col>::dtRotation(const dtQuaternion<m_type, 4> &q)
{
    Quat2RotMat(q.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 0;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 0;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetIdentity()
{
    m_elem[0] = 1;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 1;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetDiagonal(const m_type d1, const m_type d2, const m_type d3)
{
    m_elem[0] = d1;
    m_elem[4] = d2;
    m_elem[8] = d3;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetDiagonal(const m_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(m_type))
    {
    case 1:
        m_elem[0] = element[0];
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[4] = element[1];
        break;
    default:
        m_elem[0] = element[0];
        m_elem[4] = element[1];
        m_elem[8] = element[2];
        break;
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetDiagonal(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetDiagonal(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetFill(const m_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
    m_elem[4] = value;
    m_elem[5] = value;
    m_elem[6] = value;
    m_elem[7] = value;
    m_elem[8] = value;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const m_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(m_type) * m_row * m_col;

    if (matSz <= n_byte)
    {
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        m_elem[6] = element[6];
        m_elem[7] = element[7];
        m_elem[8] = element[8];
    }
    else
        memcpy(m_elem, element, n_byte);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(
    const m_type m00, const m_type m01, const m_type m02,
    const m_type m10, const m_type m11, const m_type m12,
    const m_type m20, const m_type m21, const m_type m22)
{
    m_elem[0] = m00;
    m_elem[1] = m01;
    m_elem[2] = m02;
    m_elem[3] = m10;
    m_elem[4] = m11;
    m_elem[5] = m12;
    m_elem[6] = m20;
    m_elem[7] = m21;
    m_elem[8] = m22;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const m_type angle)
{
    m_type c = std::cos(angle);
    m_type s = std::sin(angle);

    switch (order)
    {
    case 0x0: // X-axis
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = c;
        m_elem[5] = -s;
        m_elem[6] = 0;
        m_elem[7] = s;
        m_elem[8] = c;
        break;
    case 0x1: // Y-axis
        m_elem[0] = c;
        m_elem[1] = 0;
        m_elem[2] = s;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = -s;
        m_elem[7] = 0;
        m_elem[8] = c;
        break;
    case 0x2: // Z-axis
        m_elem[0] = c;
        m_elem[1] = -s;
        m_elem[2] = 0;
        m_elem[3] = s;
        m_elem[4] = c;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
        break;
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const m_type angle1, const m_type angle2)
{
    SetElement(order & 0xF, angle1); // R1
    dtRotation<m_type> R2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * R2;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3)
{
    SetElement(order & 0xF, angle1);                   // R1
    dtRotation<m_type> R2((order >> 4) & 0xF, angle2); // R2
    dtRotation<m_type> R3((order >> 8) & 0xF, angle3); // R3
    (*this) = (*this) * R2 * R3;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const dtRotation &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const dtVector3<m_type, 3> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const dtVector<3, m_type> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const uint16_t order, const m_type *e)
{
    Euler2RotMat(order, e);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const dtQuaternion<m_type, 4> &q)
{
    Quat2RotMat(q.m_elem);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const m_type *q)
{
    Quat2RotMat(q);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetElement(const m_type w, const m_type x, const m_type y, const m_type z)
{
    m_type q[4] = {w, x, y, z};
    Quat2RotMat(q);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
{
    m_type tmpElem[m_col];

    tmpElem[0] = m_elem[idxRow1 * m_col];
    m_elem[idxRow1 * m_col] = m_elem[idxRow2 * m_col];
    m_elem[idxRow2 * m_col] = tmpElem[0];

    tmpElem[1] = m_elem[idxRow1 * m_col + 1];
    m_elem[idxRow1 * m_col + 1] = m_elem[idxRow2 * m_col + 1];
    m_elem[idxRow2 * m_col + 1] = tmpElem[1];

    tmpElem[2] = m_elem[idxRow1 * m_col + 2];
    m_elem[idxRow1 * m_col + 2] = m_elem[idxRow2 * m_col + 2];
    m_elem[idxRow2 * m_col + 2] = tmpElem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
{
    m_type tmpElem[m_row];

    tmpElem[0] = m_elem[idxCol1];
    m_elem[idxCol1] = m_elem[idxCol2];
    m_elem[idxCol2] = tmpElem[0];

    tmpElem[1] = m_elem[m_col + idxCol1];
    m_elem[m_col + idxCol1] = m_elem[m_col + idxCol2];
    m_elem[m_col + idxCol2] = tmpElem[1];

    tmpElem[2] = m_elem[2 * m_col + idxCol1];
    m_elem[2 * m_col + idxCol1] = m_elem[2 * m_col + idxCol2];
    m_elem[2 * m_col + idxCol2] = tmpElem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline const m_type *const dtRotation<m_type, m_row, m_col>::GetElementsAddr() const
{
    return m_elem;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtRotation<m_type, m_row, m_col>::GetRowVec(const uint16_t idxRow) const
{
    return dtVector3<m_type, 3>(
        m_elem[idxRow * m_col],
        m_elem[idxRow * m_col + 1],
        m_elem[idxRow * m_col + 2]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtRotation<m_type, m_row, m_col>::GetColVec(const uint16_t idxCol) const
{
    return dtVector3<m_type, 3>(
        m_elem[0 * m_col + idxCol],
        m_elem[1 * m_col + idxCol],
        m_elem[2 * m_col + idxCol]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline int8_t dtRotation<m_type, m_row, m_col>::GetRowVec(const uint16_t idxRow, dtVector3<m_type, 3> &v) const
{
    v.m_elem[0] = m_elem[idxRow * m_col];
    v.m_elem[1] = m_elem[idxRow * m_col + 1];
    v.m_elem[2] = m_elem[idxRow * m_col + 2];

    return 0;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline int8_t dtRotation<m_type, m_row, m_col>::GetColVec(const uint16_t idxCol, dtVector3<m_type, 3> &v) const
{
    v.m_elem[0] = m_elem[0 * m_col + idxCol];
    v.m_elem[1] = m_elem[1 * m_col + idxCol];
    v.m_elem[2] = m_elem[2 * m_col + idxCol];

    return 0;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtRotation<m_type, m_row, m_col>::GetEulerAngles(uint16_t order) const
{
    /* Tait?Bryan angles */
    m_type vec[3];

    // 0:x, 1:y, 2:z
    // order = 0x012 -> zyx, 0x210 -> xyz, 0x102 -> zxy, inverse order!
    uint16_t o1 = order & 0xF;
    uint16_t o2 = (order >> 4) & 0xF;
    uint16_t o3 = (order >> 8) & 0xF;
    int sign = ((o1 + 1) == o2) ? 1 : -1;

    vec[1] = std::asin(sign * m_elem[o1 * m_col + o3]);

    if ((1 - std::fabs(m_elem[o1 * m_col + o3])) <= std::numeric_limits<m_type>::epsilon())
    {
        // case s2 = +-1, c2 = 0 s3 = 0 c3 = +-1
        vec[0] = std::atan2(sign * m_elem[o3 * m_col + o2], m_elem[o2 * m_col + o2]);
        vec[2] = 0;
    }
    else
    {
        // case s2 != 1 or s2 != -1
        vec[0] = std::atan2(-sign * m_elem[o2 * m_col + o3], m_elem[o3 * m_col + o3]);
        vec[2] = std::atan2(-sign * m_elem[o1 * m_col + o2], m_elem[o1 * m_col + o1]);
    }

    return dtVector3<m_type, 3>(vec);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::Transpose() const
{
    return dtRotation(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::log() const
{
    /* Logarithm of rotation matrix */
    // That is Axis-angle representation
    // log:SO(3) -> so(3); R -> log(R) = [u*phi]x
    // phi = arccos((trace(R)-1)/2)
    // u = (R-RT)V / 2sin(phi)
    // ()V is the inverse of []x
    // if phi is 0, singular

    m_type phi = std::acos((m_elem[0] + m_elem[4] + m_elem[8] - 1) * (m_type)(0.5));
    m_type alpha;

    if (std::abs(phi) > std::numeric_limits<m_type>::epsilon())
        alpha = phi / (2 * std::sin(phi));
    else
        alpha = 0;

    return dtRotation(
        0, (m_elem[1] - m_elem[3]) * alpha, (m_elem[2] - m_elem[6]) * alpha,
        (m_elem[3] - m_elem[1]) * alpha, 0, (m_elem[5] - m_elem[7]) * alpha,
        (m_elem[6] - m_elem[2]) * alpha, (m_elem[7] - m_elem[5]) * alpha, 0);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtRotation<m_type, m_row, m_col>::Log() const
{
    /* Logarithm of rotation matrix */
    // That is Axis-angle representation
    // Log:SO(3) -> R3; R -> Log(R) = u*phi
    // phi = arccos((trace(R)-1)/2)
    // u = (R-RT)V / 2sin(phi)
    // ()V is the inverse of []x, matrix([]x) -> vector(()V)
    // if phi is 0, singular

    m_type phi = std::acos((m_elem[0] + m_elem[4] + m_elem[8] - 1) * (m_type)(0.5));
    m_type alpha;

    if (std::abs(phi) > std::numeric_limits<m_type>::epsilon())
        alpha = phi / (2 * std::sin(phi));
    else
        alpha = 0;

    return dtVector3<m_type, 3>(
        (m_elem[7] - m_elem[5]) * alpha,
        (m_elem[2] - m_elem[6]) * alpha,
        (m_elem[3] - m_elem[1]) * alpha);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::ode(m_type wx, m_type wy, m_type wz) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return dtRotation(
        m_elem[1] * wz - m_elem[2] * wy, -m_elem[0] * wz + m_elem[2] * wx, m_elem[0] * wy - m_elem[1] * wx,
        m_elem[4] * wz - m_elem[5] * wy, -m_elem[3] * wz + m_elem[5] * wx, m_elem[3] * wy - m_elem[4] * wx,
        m_elem[7] * wz - m_elem[8] * wy, -m_elem[6] * wz + m_elem[8] * wx, m_elem[6] * wy - m_elem[7] * wx);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::ode(m_type *w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return dtRotation(
        m_elem[1] * w[2] - m_elem[2] * w[1], -m_elem[0] * w[2] + m_elem[2] * w[0], m_elem[0] * w[1] - m_elem[1] * w[0],
        m_elem[4] * w[2] - m_elem[5] * w[1], -m_elem[3] * w[2] + m_elem[5] * w[0], m_elem[3] * w[1] - m_elem[4] * w[0],
        m_elem[7] * w[2] - m_elem[8] * w[1], -m_elem[6] * w[2] + m_elem[8] * w[0], m_elem[6] * w[1] - m_elem[7] * w[0]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::ode(dtVector3<m_type, 3> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return dtRotation(
        m_elem[1] * w.m_elem[2] - m_elem[2] * w.m_elem[1], -m_elem[0] * w.m_elem[2] + m_elem[2] * w.m_elem[0], m_elem[0] * w.m_elem[1] - m_elem[1] * w.m_elem[0],
        m_elem[4] * w.m_elem[2] - m_elem[5] * w.m_elem[1], -m_elem[3] * w.m_elem[2] + m_elem[5] * w.m_elem[0], m_elem[3] * w.m_elem[1] - m_elem[4] * w.m_elem[0],
        m_elem[7] * w.m_elem[2] - m_elem[8] * w.m_elem[1], -m_elem[6] * w.m_elem[2] + m_elem[8] * w.m_elem[0], m_elem[6] * w.m_elem[1] - m_elem[7] * w.m_elem[0]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::ode(dtVector<3, m_type> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return dtRotation(
        m_elem[1] * w.m_elem[2] - m_elem[2] * w.m_elem[1], -m_elem[0] * w.m_elem[2] + m_elem[2] * w.m_elem[0], m_elem[0] * w.m_elem[1] - m_elem[1] * w.m_elem[0],
        m_elem[4] * w.m_elem[2] - m_elem[5] * w.m_elem[1], -m_elem[3] * w.m_elem[2] + m_elem[5] * w.m_elem[0], m_elem[3] * w.m_elem[1] - m_elem[4] * w.m_elem[0],
        m_elem[7] * w.m_elem[2] - m_elem[8] * w.m_elem[1], -m_elem[6] * w.m_elem[2] + m_elem[8] * w.m_elem[0], m_elem[6] * w.m_elem[1] - m_elem[7] * w.m_elem[0]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::Inv() const
{
    return dtRotation(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

/* Assignment operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> &dtRotation<m_type, m_row, m_col>::operator=(const dtRotation &m)
{
    // memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtCommaInit<m_row * m_col, m_type> dtRotation<m_type, m_row, m_col>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row * m_col, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator-() const
{
    return dtRotation(
        -m_elem[0], -m_elem[1], -m_elem[2],
        -m_elem[3], -m_elem[4], -m_elem[5],
        -m_elem[6], -m_elem[7], -m_elem[8]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator+(const dtRotation &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator-(const dtRotation &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator+(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator-(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator+(const dtMatrix<m_row, m_col, m_type> &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator-(const dtMatrix<m_row, m_col, m_type> &m) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator*(const m_type s) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) *= s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator/(const m_type s) const
{
    return dtMatrix3<m_type, m_row, m_col>(*this) /= s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtRotation<m_type, m_row, m_col>::operator*(const dtMatrix<m_col, col, m_type> &m) const
{
    m_type mat[m_row * col];

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        for (uint16_t icol = 0; icol < col; ++icol)
        {
            mat[irow * col + icol] = m_elem[irow * m_col] * m.m_elem[icol];
            mat[irow * col + icol] += m_elem[irow * m_col + 1] * m.m_elem[col + icol];
            mat[irow * col + icol] += m_elem[irow * m_col + 2] * m.m_elem[2 * col + icol];
        }
    }

    return dtMatrix<m_row, col, m_type>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator*(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] * m.m_elem[0] + m_elem[1] * m.m_elem[3] + m_elem[2] * m.m_elem[6];
    mat[1] = m_elem[0] * m.m_elem[1] + m_elem[1] * m.m_elem[4] + m_elem[2] * m.m_elem[7];
    mat[2] = m_elem[0] * m.m_elem[2] + m_elem[1] * m.m_elem[5] + m_elem[2] * m.m_elem[8];

    mat[3] = m_elem[3] * m.m_elem[0] + m_elem[4] * m.m_elem[3] + m_elem[5] * m.m_elem[6];
    mat[4] = m_elem[3] * m.m_elem[1] + m_elem[4] * m.m_elem[4] + m_elem[5] * m.m_elem[7];
    mat[5] = m_elem[3] * m.m_elem[2] + m_elem[4] * m.m_elem[5] + m_elem[5] * m.m_elem[8];

    mat[6] = m_elem[6] * m.m_elem[0] + m_elem[7] * m.m_elem[3] + m_elem[8] * m.m_elem[6];
    mat[7] = m_elem[6] * m.m_elem[1] + m_elem[7] * m.m_elem[4] + m_elem[8] * m.m_elem[7];
    mat[8] = m_elem[6] * m.m_elem[2] + m_elem[7] * m.m_elem[5] + m_elem[8] * m.m_elem[8];

    return dtMatrix3<m_type, m_row, m_col>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator*(const dtRotation &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] * m.m_elem[0] + m_elem[1] * m.m_elem[3] + m_elem[2] * m.m_elem[6];
    mat[1] = m_elem[0] * m.m_elem[1] + m_elem[1] * m.m_elem[4] + m_elem[2] * m.m_elem[7];
    mat[2] = m_elem[0] * m.m_elem[2] + m_elem[1] * m.m_elem[5] + m_elem[2] * m.m_elem[8];

    mat[3] = m_elem[3] * m.m_elem[0] + m_elem[4] * m.m_elem[3] + m_elem[5] * m.m_elem[6];
    mat[4] = m_elem[3] * m.m_elem[1] + m_elem[4] * m.m_elem[4] + m_elem[5] * m.m_elem[7];
    mat[5] = m_elem[3] * m.m_elem[2] + m_elem[4] * m.m_elem[5] + m_elem[5] * m.m_elem[8];

    mat[6] = m_elem[6] * m.m_elem[0] + m_elem[7] * m.m_elem[3] + m_elem[8] * m.m_elem[6];
    mat[7] = m_elem[6] * m.m_elem[1] + m_elem[7] * m.m_elem[4] + m_elem[8] * m.m_elem[7];
    mat[8] = m_elem[6] * m.m_elem[2] + m_elem[7] * m.m_elem[5] + m_elem[8] * m.m_elem[8];

    return dtRotation(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector<m_row, m_type> dtRotation<m_type, m_row, m_col>::operator*(const dtVector<m_col, m_type> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return dtVector<m_row, m_type>(vec);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, m_row> dtRotation<m_type, m_row, m_col>::operator*(const dtVector3<m_type, m_col> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return dtVector3<m_type, m_row>(vec);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator&(const dtVector<m_col, m_type> &v) const
{ // RotMat * [v]x, []x is skew-symmetric matrix
    return dtRotation(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtRotation<m_type, m_row, m_col> dtRotation<m_type, m_row, m_col>::operator&(const dtVector3<m_type, m_col> &v) const
{ // RotMat * [v]x, []x is skew-symmetric matrix
    return dtRotation(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator==(const dtRotation &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator!=(const dtRotation &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator==(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator!=(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator==(const dtMatrix<m_row, m_col, m_type> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtRotation<m_type, m_row, m_col>::operator!=(const dtMatrix<m_row, m_col, m_type> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance)
        return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance)
        return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance)
        return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
        {
            Serial.printf("%7.3f ", (m_type)(m_elem[irow * m_col + icol]));
        }
        Serial.write('\n');
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
        {
            printf("%7.3f ", (m_type)(m_elem[irow * m_col + icol]));
        }
        printf("\n");
    }
    printf("%c", endChar);
#endif
}

//-- Private Member Function ------------------------------------------------//
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::Euler2RotMat(const uint16_t order, const m_type *e)
{
    m_type s1 = std::sin(e[0]);
    m_type c1 = std::cos(e[0]);
    m_type s2 = std::sin(e[1]);
    m_type c2 = std::cos(e[1]);
    m_type s3 = std::sin(e[2]);
    m_type c3 = std::cos(e[2]);

    /* Only Tait?Bryan angles */
    switch (order)
    {
    case 0x120: // xzy
        m_elem[0] = c2 * c3;
        m_elem[1] = -s2;
        m_elem[2] = c2 * s3;
        m_elem[3] = s1 * s3 + c1 * c3 * s2;
        m_elem[4] = c1 * c2;
        m_elem[5] = c1 * s2 * s3 - c3 * s1;
        m_elem[6] = c3 * s1 * s2 - c1 * s3;
        m_elem[7] = c2 * s1;
        m_elem[8] = c1 * c3 + s1 * s2 * s3;
        break;
    case 0x210: // xyz
        m_elem[0] = c2 * c3;
        m_elem[1] = -c2 * s3;
        m_elem[2] = s2;
        m_elem[3] = c1 * s3 + c3 * s1 * s2;
        m_elem[4] = c1 * c3 - s1 * s2 * s3;
        m_elem[5] = -c2 * s1;
        m_elem[6] = s1 * s3 - c1 * c3 * s2;
        m_elem[7] = c3 * s1 + c1 * s2 * s3;
        m_elem[8] = c1 * c2;
        break;
    case 0x201: // yxz
        m_elem[0] = c1 * c3 + s1 * s2 * s3;
        m_elem[1] = c3 * s1 * s2 - c1 * s3;
        m_elem[2] = c2 * s1;
        m_elem[3] = c2 * s3;
        m_elem[4] = c2 * c3;
        m_elem[5] = -s2;
        m_elem[6] = c1 * s2 * s3 - c3 * s1;
        m_elem[7] = c1 * c3 * s2 + s1 * s3;
        m_elem[8] = c1 * c2;
        break;
    case 0x021: // yzx
        m_elem[0] = c1 * c2;
        m_elem[1] = s1 * s3 - c1 * c3 * s2;
        m_elem[2] = c3 * s1 + c1 * s2 * s3;
        m_elem[3] = s2;
        m_elem[4] = c2 * c3;
        m_elem[5] = -c2 * s3;
        m_elem[6] = -c2 * s1;
        m_elem[7] = c1 * s3 + c3 * s1 * s2;
        m_elem[8] = c1 * c3 - s1 * s2 * s3;
        break;
    case 0x012: // zyx
        m_elem[0] = c1 * c2;
        m_elem[1] = c1 * s2 * s3 - c3 * s1;
        m_elem[2] = s1 * s3 + c1 * c3 * s2;
        m_elem[3] = c2 * s1;
        m_elem[4] = c1 * c3 + s1 * s2 * s3;
        m_elem[5] = c3 * s1 * s2 - c1 * s3;
        m_elem[6] = -s2;
        m_elem[7] = c2 * s3;
        m_elem[8] = c2 * c3;
        break;
    case 0x102: // zxy
        m_elem[0] = c1 * c3 - s1 * s2 * s3;
        m_elem[1] = -c2 * s1;
        m_elem[2] = c1 * s3 + c3 * s1 * s2;
        m_elem[3] = c3 * s1 + c1 * s2 * s3;
        m_elem[4] = c1 * c2;
        m_elem[5] = s1 * s3 - c1 * c3 * s2;
        m_elem[6] = -c2 * s3;
        m_elem[7] = s2;
        m_elem[8] = c2 * c3;
        break;
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtRotation<m_type, m_row, m_col>::Quat2RotMat(const m_type *q)
{
    // m_elem[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    m_elem[0] = 1 - 2 * (q[2] * q[2] + q[3] * q[3]);
    m_elem[1] = 2 * (q[1] * q[2] - q[0] * q[3]);
    m_elem[2] = 2 * (q[1] * q[3] + q[0] * q[2]);

    m_elem[3] = 2 * (q[1] * q[2] + q[0] * q[3]);
    // m_elem[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    m_elem[4] = 1 - 2 * (q[1] * q[1] + q[3] * q[3]);
    m_elem[5] = 2 * (q[2] * q[3] - q[0] * q[1]);

    m_elem[6] = 2 * (q[1] * q[3] - q[0] * q[2]);
    m_elem[7] = 2 * (q[2] * q[3] + q[0] * q[1]);
    // m_elem[8] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
    m_elem[8] = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <typename type, uint16_t row, uint16_t col>
inline dtMatrix3<type, row, col> operator*(const type s, const dtRotation<type, row, col> &m)
{
    return dtMatrix3<type, row, col>(m) *= s;
}

typedef dtRotation<> dtRotMat;

} // namespace dtMath

#endif // DTMATH_DTROTATION_TPP_