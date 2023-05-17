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

#ifndef DTMATH_DTVECTOR6_TPP_
#define DTMATH_DTVECTOR6_TPP_

#include "dtVector6.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 0;
    m_elem[5] = 0;
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const m_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
    m_elem[4] = element[4];
    m_elem[5] = element[5];
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const m_type *element, const size_t n_byte)
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
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = 0;
        m_elem[5] = 0;
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = 0;
        m_elem[4] = 0;
        m_elem[5] = 0;
        break;
    case 4:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = 0;
        m_elem[5] = 0;
        break;
    case 5:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = 0;
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz)
{
    m_elem[0] = px;
    m_elem[1] = py;
    m_elem[2] = pz;
    m_elem[3] = ox;
    m_elem[4] = oy;
    m_elem[5] = oz;
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const dtVector6<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const dtVector3<m_type, 3> &p, const dtVector3<m_type, 3> &o)
{
    m_elem[0] = p.m_elem[0];
    m_elem[1] = p.m_elem[1];
    m_elem[2] = p.m_elem[2];
    m_elem[3] = o.m_elem[0];
    m_elem[4] = o.m_elem[1];
    m_elem[5] = o.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row>::dtVector6(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 0;
    m_elem[5] = 0;
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetFill(const m_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
    m_elem[4] = value;
    m_elem[5] = value;
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const m_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(m_type))
    {
    case 1:
        m_elem[0] = element[0];
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        break;
    case 4:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        break;
    case 5:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz)
{
    m_elem[0] = px;
    m_elem[1] = py;
    m_elem[2] = pz;
    m_elem[3] = ox;
    m_elem[4] = oy;
    m_elem[5] = oz;
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const dtVector6 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const dtVector3<m_type, 3> &p, const dtVector3<m_type, 3> &o)
{
    m_elem[0] = p.m_elem[0];
    m_elem[1] = p.m_elem[1];
    m_elem[2] = p.m_elem[2];
    m_elem[3] = o.m_elem[0];
    m_elem[4] = o.m_elem[1];
    m_elem[5] = o.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetElement(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector<row, m_type> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row)
        rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    case 4:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    case 5:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        m_elem[idxRow + 4] = v.m_elem[4];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        m_elem[4] = v.m_elem[4];
        m_elem[5] = v.m_elem[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const m_type *v, const size_t n_byte)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    uint16_t row = n_byte / sizeof(m_type);
    if (rowSz > row)
        rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v[0];
        break;
    case 2:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        break;
    case 3:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        m_elem[idxRow + 2] = v[2];
        break;
    case 4:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        m_elem[idxRow + 2] = v[2];
        m_elem[idxRow + 3] = v[3];
        break;
    case 5:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        m_elem[idxRow + 2] = v[2];
        m_elem[idxRow + 3] = v[3];
        m_elem[idxRow + 4] = v[4];
        break;
    default:
        m_elem[0] = v[0];
        m_elem[1] = v[1];
        m_elem[2] = v[2];
        m_elem[3] = v[3];
        m_elem[4] = v[4];
        m_elem[5] = v[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector3<m_type, 3> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > 3)
        rowSz = 3;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector4<m_type, 4> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > 4)
        rowSz = 4;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    case 4:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector6<m_type, 6> &v)
{
    if (idxRow >= m_row)
        return;

    switch (idxRow)
    {
    case 5:
        m_elem[5] = v.m_elem[0];
        break;
    case 4:
        m_elem[4] = v.m_elem[0];
        m_elem[5] = v.m_elem[1];
        break;
    case 3:
        m_elem[3] = v.m_elem[0];
        m_elem[4] = v.m_elem[1];
        m_elem[5] = v.m_elem[2];
        break;
    case 2:
        m_elem[2] = v.m_elem[0];
        m_elem[3] = v.m_elem[1];
        m_elem[4] = v.m_elem[2];
        m_elem[5] = v.m_elem[3];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        m_elem[3] = v.m_elem[2];
        m_elem[4] = v.m_elem[3];
        m_elem[5] = v.m_elem[4];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        m_elem[4] = v.m_elem[4];
        m_elem[5] = v.m_elem[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline void dtVector6<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtMatrix<row, 1, m_type> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row)
        rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    case 4:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    case 5:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        m_elem[idxRow + 4] = v.m_elem[4];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        m_elem[4] = v.m_elem[4];
        m_elem[5] = v.m_elem[5];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetSwap(const uint16_t i, const uint16_t j)
{
    m_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::SetNormalize()
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] +
        m_elem[4] * m_elem[4] +
        m_elem[5] * m_elem[5]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    m_elem[0] /= norm;
    m_elem[1] /= norm;
    m_elem[2] /= norm;
    m_elem[3] /= norm;
    m_elem[4] /= norm;
    m_elem[5] /= norm;
}

template <typename m_type, uint16_t m_row>
inline const m_type *const dtVector6<m_type, m_row>::GetElementsAddr() const
{
    return m_elem;
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline dtVector<row, m_type> dtVector6<m_type, m_row>::GetBlock(const uint16_t idx)
{
    m_type elem[row] = {
        0,
    };
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row)
        return dtVector<row, m_type>(elem);
    if (rowSize > row)
        rowSize = row;

    memcpy(elem, &m_elem[idx], sizeof(m_type) * rowSize);

    return dtVector<row, m_type>(elem);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtVector6<m_type, m_row>::GetBlockVec3(const uint16_t idx)
{
    m_type elem[3] = {
        0,
    };

    if (idx >= m_row)
        return dtVector3<m_type, 3>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
    };

    return dtVector3<m_type, 3>(elem);
}

template <typename m_type, uint16_t m_row>
inline dtVector4<m_type, 4> dtVector6<m_type, m_row>::GetBlockVec4(const uint16_t idx)
{
    m_type elem[4] = {
        0,
    };

    if (idx >= m_row)
        return dtVector4<m_type, 4>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    case 3:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
    };

    return dtVector4<m_type, 4>(elem);
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline int8_t dtVector6<m_type, m_row>::GetBlock(const uint16_t idx, dtVector<row, m_type> &v)
{
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row)
        return -1;
    if (rowSize > row)
        rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(m_type) * rowSize);

    return 0;
}

template <typename m_type, uint16_t m_row>
inline int8_t dtVector6<m_type, m_row>::GetBlockVec3(const uint16_t idx, dtVector3<m_type, 3> &v)
{
    if (idx >= m_row)
        return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
    };

    return 0;
}

template <typename m_type, uint16_t m_row>
inline int8_t dtVector6<m_type, m_row>::GetBlockVec4(const uint16_t idx, dtVector4<m_type, 4> &v)
{
    if (idx >= m_row)
        return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
    };

    return 0;
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtVector6<m_type, m_row>::GetPos() const
{
    return dtVector3<m_type, 3>(m_elem[0], m_elem[1], m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtVector6<m_type, m_row>::GetOri() const
{
    return dtVector3<m_type, 3>(m_elem[3], m_elem[4], m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::GetNorm() const
{
    return std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] +
        m_elem[4] * m_elem[4] +
        m_elem[5] * m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::GetSqNorm() const
{
    return (
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] +
        m_elem[4] * m_elem[4] +
        m_elem[5] * m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::GetSum() const
{
    return (
        m_elem[0] +
        m_elem[1] +
        m_elem[2] +
        m_elem[3] +
        m_elem[4] +
        m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::GetNormalized() const
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] +
        m_elem[4] * m_elem[4] +
        m_elem[5] * m_elem[5]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    return dtVector6(
        m_elem[0] / norm,
        m_elem[1] / norm,
        m_elem[2] / norm,
        m_elem[3] / norm,
        m_elem[4] / norm,
        m_elem[5] / norm);
}

template <typename m_type, uint16_t m_row>
inline dtMatrix<1, m_row, m_type> dtVector6<m_type, m_row>::Transpose() const
{
    return dtMatrix<1, m_row, m_type>(m_elem);
}

/* Assignment operators */
template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator=(const dtVector6 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator+=(const dtVector6 &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator-=(const dtVector6 &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator*=(const dtVector6 &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator/=(const dtVector6 &v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[4] /= den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[5] /= den;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator+=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator-=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator*=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator/=(const dtVector<m_row, m_type> &v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[4] /= den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[5] /= den;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator+=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator-=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator*=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator/=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[3] /= den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[4] /= den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    m_elem[5] /= den;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator=(const m_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;
    m_elem[3] = s;
    m_elem[4] = s;
    m_elem[5] = s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator+=(const m_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;
    m_elem[3] += s;
    m_elem[4] += s;
    m_elem[5] += s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator-=(const m_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;
    m_elem[3] -= s;
    m_elem[4] -= s;
    m_elem[5] -= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator*=(const m_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;
    m_elem[4] *= s;
    m_elem[5] *= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> &dtVector6<m_type, m_row>::operator/=(const m_type s)
{
    m_type den = s;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }

    m_elem[0] /= den;
    m_elem[1] /= den;
    m_elem[2] /= den;
    m_elem[3] /= den;
    m_elem[4] /= den;
    m_elem[5] /= den;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtCommaInit<m_row, m_type> dtVector6<m_type, m_row>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator-() const
{
    return dtVector6(
        -m_elem[0],
        -m_elem[1],
        -m_elem[2],
        -m_elem[3],
        -m_elem[4],
        -m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator+(const dtVector6 &v) const
{
    return dtVector6(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3],
        m_elem[4] + v.m_elem[4],
        m_elem[5] + v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator-(const dtVector6 &v) const
{
    return dtVector6(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3],
        m_elem[4] - v.m_elem[4],
        m_elem[5] - v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator*(const dtVector6 &v) const
{
    return dtVector6(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3],
        m_elem[4] * v.m_elem[4],
        m_elem[5] * v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator/(const dtVector6 &v) const
{
    m_type den[6];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];
    den[4] = v.m_elem[4];
    den[5] = v.m_elem[5];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0)
            den[0] = -std::numeric_limits<m_type>::epsilon();
        else
            den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0)
            den[1] = -std::numeric_limits<m_type>::epsilon();
        else
            den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0)
            den[2] = -std::numeric_limits<m_type>::epsilon();
        else
            den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0)
            den[3] = -std::numeric_limits<m_type>::epsilon();
        else
            den[3] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[4]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[4] < 0)
            den[4] = -std::numeric_limits<m_type>::epsilon();
        else
            den[4] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[5]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[5] < 0)
            den[5] = -std::numeric_limits<m_type>::epsilon();
        else
            den[5] = std::numeric_limits<m_type>::epsilon();
    }

    return dtVector6(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3],
        m_elem[4] / den[4],
        m_elem[5] / den[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator+(const dtVector<m_row, m_type> &v) const
{
    return dtVector6(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3],
        m_elem[4] + v.m_elem[4],
        m_elem[5] + v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator-(const dtVector<m_row, m_type> &v) const
{
    return dtVector6(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3],
        m_elem[4] - v.m_elem[4],
        m_elem[5] - v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator*(const dtVector<m_row, m_type> &v) const
{
    return dtVector6(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3],
        m_elem[4] * v.m_elem[4],
        m_elem[5] * v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator/(const dtVector<m_row, m_type> &v) const
{
    m_type den[6];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];
    den[4] = v.m_elem[4];
    den[5] = v.m_elem[5];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0)
            den[0] = -std::numeric_limits<m_type>::epsilon();
        else
            den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0)
            den[1] = -std::numeric_limits<m_type>::epsilon();
        else
            den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0)
            den[2] = -std::numeric_limits<m_type>::epsilon();
        else
            den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0)
            den[3] = -std::numeric_limits<m_type>::epsilon();
        else
            den[3] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[4]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[4] < 0)
            den[4] = -std::numeric_limits<m_type>::epsilon();
        else
            den[4] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[5]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[5] < 0)
            den[5] = -std::numeric_limits<m_type>::epsilon();
        else
            den[5] = std::numeric_limits<m_type>::epsilon();
    }

    return dtVector6(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3],
        m_elem[4] / den[4],
        m_elem[5] / den[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator+(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector6(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3],
        m_elem[4] + v.m_elem[4],
        m_elem[5] + v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator-(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector6(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3],
        m_elem[4] - v.m_elem[4],
        m_elem[5] - v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator*(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector6(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3],
        m_elem[4] * v.m_elem[4],
        m_elem[5] * v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator/(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type den[6];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];
    den[4] = v.m_elem[4];
    den[5] = v.m_elem[5];

    if (std::abs(den[0]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[0] < 0)
            den[0] = -std::numeric_limits<m_type>::epsilon();
        else
            den[0] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[1] < 0)
            den[1] = -std::numeric_limits<m_type>::epsilon();
        else
            den[1] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[2] < 0)
            den[2] = -std::numeric_limits<m_type>::epsilon();
        else
            den[2] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[3] < 0)
            den[3] = -std::numeric_limits<m_type>::epsilon();
        else
            den[3] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[4]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[4] < 0)
            den[4] = -std::numeric_limits<m_type>::epsilon();
        else
            den[4] = std::numeric_limits<m_type>::epsilon();
    }
    if (std::abs(den[5]) < std::numeric_limits<m_type>::epsilon())
    {
        if (den[5] < 0)
            den[5] = -std::numeric_limits<m_type>::epsilon();
        else
            den[5] = std::numeric_limits<m_type>::epsilon();
    }

    return dtVector6(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3],
        m_elem[4] / den[4],
        m_elem[5] / den[5]);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator+(const m_type s) const
{
    return dtVector6(
        m_elem[0] + s,
        m_elem[1] + s,
        m_elem[2] + s,
        m_elem[3] + s,
        m_elem[4] + s,
        m_elem[5] + s);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator-(const m_type s) const
{
    return dtVector6(
        m_elem[0] - s,
        m_elem[1] - s,
        m_elem[2] - s,
        m_elem[3] - s,
        m_elem[4] - s,
        m_elem[5] - s);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator*(const m_type s) const
{
    return dtVector6(
        m_elem[0] * s,
        m_elem[1] * s,
        m_elem[2] * s,
        m_elem[3] * s,
        m_elem[4] * s,
        m_elem[5] * s);
}

template <typename m_type, uint16_t m_row>
inline dtVector6<m_type, m_row> dtVector6<m_type, m_row>::operator/(const m_type s) const
{
    m_type den = s;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }

    return dtVector6(
        m_elem[0] / den,
        m_elem[1] / den,
        m_elem[2] / den,
        m_elem[3] / den,
        m_elem[4] / den,
        m_elem[5] / den);
}

template <typename m_type, uint16_t m_row>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtVector6<m_type, m_row>::operator*(const dtMatrix<1, col, m_type> &m) const
{
    m_type mat[m_row * col];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
            mat[irow * col + icol + 1] = m_elem[irow] * m.m_elem[icol + 1];
            mat[irow * col + icol + 2] = m_elem[irow] * m.m_elem[icol + 2];
            mat[irow * col + icol + 3] = m_elem[irow] * m.m_elem[icol + 3];
        }

        for (cnt = col % 4u; cnt > 0u; cnt--, icol++)
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
    }

    return dtMatrix<m_row, col, m_type>(mat);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::dot(const dtVector6 &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::dot(const dtVector<m_row, m_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector6<m_type, m_row>::dot(const dtMatrix<m_row, 1, m_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator==(const dtVector6 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator!=(const dtVector6 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator==(const dtVector<m_row, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator!=(const dtVector<m_row, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator==(const dtMatrix<m_row, 1, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector6<m_type, m_row>::operator!=(const dtMatrix<m_row, 1, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline void dtVector6<m_type, m_row>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        Serial.printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <typename type, uint16_t row>
inline dtVector6<type, row> operator+(const type s, const dtVector6<type, row> &v)
{
    return dtVector6<type, row>(
        v.m_elem[0] + s,
        v.m_elem[1] + s,
        v.m_elem[2] + s,
        v.m_elem[3] + s,
        v.m_elem[4] + s,
        v.m_elem[5] + s);
}

// scalar - vector
template <typename type, uint16_t row>
inline dtVector6<type, row> operator-(const type s, const dtVector6<type, row> &v)
{
    return dtVector6<type, row>(
        s - v.m_elem[0],
        s - v.m_elem[1],
        s - v.m_elem[2],
        s - v.m_elem[3],
        s - v.m_elem[4],
        s - v.m_elem[5]);
}

// scalar * vector
template <typename type, uint16_t row>
inline dtVector6<type, row> operator*(const type s, const dtVector6<type, row> &v)
{
    return dtVector6<type, row>(
        v.m_elem[0] * s,
        v.m_elem[1] * s,
        v.m_elem[2] * s,
        v.m_elem[3] * s,
        v.m_elem[4] * s,
        v.m_elem[5] * s);
}

// scalar / vector
template <typename type, uint16_t row>
inline dtVector6<type, row> operator/(const type s, const dtVector6<type, row> &v)
{
    type den[6];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];
    den[4] = v.m_elem[4];
    den[5] = v.m_elem[5];

    if (std::abs(den[0]) < std::numeric_limits<type>::epsilon())
    {
        if (den[0] < 0)
            den[0] = -std::numeric_limits<type>::epsilon();
        else
            den[0] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<type>::epsilon())
    {
        if (den[1] < 0)
            den[1] = -std::numeric_limits<type>::epsilon();
        else
            den[1] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<type>::epsilon())
    {
        if (den[2] < 0)
            den[2] = -std::numeric_limits<type>::epsilon();
        else
            den[2] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<type>::epsilon())
    {
        if (den[3] < 0)
            den[3] = -std::numeric_limits<type>::epsilon();
        else
            den[3] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[4]) < std::numeric_limits<type>::epsilon())
    {
        if (den[4] < 0)
            den[4] = -std::numeric_limits<type>::epsilon();
        else
            den[4] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[5]) < std::numeric_limits<type>::epsilon())
    {
        if (den[5] < 0)
            den[5] = -std::numeric_limits<type>::epsilon();
        else
            den[5] = std::numeric_limits<type>::epsilon();
    }

    return dtVector6<type, row>(
        s / den[0],
        s / den[1],
        s / den[2],
        s / den[3],
        s / den[4],
        s / den[5]);
}

typedef dtVector6<> dtVec6;

} // namespace dtMath

#endif // DTMATH_DTVECTOR6_TPP_
