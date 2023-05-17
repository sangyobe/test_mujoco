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

#ifndef DTMATH_DTVECTOR3_TPP_
#define DTMATH_DTVECTOR3_TPP_

#include "dtVector3.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const m_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const m_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(m_type))
    {
    case 1:
        m_elem[0] = element[0];
        m_elem[1] = 0;
        m_elem[2] = 0;
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = 0;
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const m_type i, const m_type j, const m_type k)
{
    m_elem[0] = i;
    m_elem[1] = j;
    m_elem[2] = k;
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const dtVector3 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row>::dtVector3(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetFill(const m_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetElement(const m_type *element, const size_t n_byte)
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
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetElement(const m_type i, const m_type j, const m_type k)
{
    m_elem[0] = i;
    m_elem[1] = j;
    m_elem[2] = k;
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetElement(const dtVector3 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetElement(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetElement(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector<row, m_type> &v)
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
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const m_type *v, const size_t n_byte)
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
    default:
        m_elem[0] = v[0];
        m_elem[1] = v[1];
        m_elem[2] = v[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector3<m_type, 3> &v)
{
    if (idxRow >= m_row)
        return;

    switch (idxRow)
    {
    case 2:
        m_elem[2] = v.m_elem[0];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector4<m_type, 4> &v)
{
    if (idxRow >= m_row)
        return;

    switch (idxRow)
    {
    case 2:
        m_elem[2] = v.m_elem[0];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtVector6<m_type, 6> &v)
{
    if (idxRow >= m_row)
        return;

    switch (idxRow)
    {
    case 2:
        m_elem[2] = v.m_elem[0];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline void dtVector3<m_type, m_row>::SetBlock(const uint16_t idxRow, const dtMatrix<row, 1, m_type> &v)
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
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetSwap(const uint16_t i, const uint16_t j)
{
    m_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::SetNormalize()
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    m_elem[0] /= norm;
    m_elem[1] /= norm;
    m_elem[2] /= norm;
}

template <typename m_type, uint16_t m_row>
inline const m_type *const dtVector3<m_type, m_row>::GetElementsAddr() const
{
    return m_elem;
}

template <typename m_type, uint16_t m_row>
template <uint16_t row>
inline dtVector<row, m_type> dtVector3<m_type, m_row>::GetBlock(const uint16_t idx)
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
template <uint16_t row>
inline int8_t dtVector3<m_type, m_row>::GetBlock(const uint16_t idx, dtVector<row, m_type> &v)
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
inline m_type dtVector3<m_type, m_row>::GetNorm() const
{
    return std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector3<m_type, m_row>::GetSqNorm() const
{
    return (
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector3<m_type, m_row>::GetSum() const
{
    return (
        m_elem[0] +
        m_elem[1] +
        m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::GetNormalized() const
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    return dtVector3<m_type, m_row>(
        m_elem[0] / norm,
        m_elem[1] / norm,
        m_elem[2] / norm);
}

template <typename m_type, uint16_t m_row>
inline dtMatrix3<m_type, m_row, m_row> dtVector3<m_type, m_row>::GetSkew() const
{
    return dtMatrix3<m_type, m_row, m_row>(
        0, -m_elem[2], m_elem[1],
        m_elem[2], 0, -m_elem[0],
        -m_elem[1], m_elem[0], 0);
}

template <typename m_type, uint16_t m_row>
inline dtMatrix<1, m_row, m_type> dtVector3<m_type, m_row>::Transpose() const
{
    return dtMatrix<1, m_row, m_type>(m_elem);
}

/* Assignment operators */
template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator=(const dtVector3 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator+=(const dtVector3 &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator-=(const dtVector3 &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator*=(const dtVector3 &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator/=(const dtVector3 &v)
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

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator+=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator-=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator*=(const dtVector<m_row, m_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator/=(const dtVector<m_row, m_type> &v)
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

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator+=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator-=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator*=(const dtMatrix<m_row, 1, m_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator/=(const dtMatrix<m_row, 1, m_type> &v)
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

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator=(const m_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator+=(const m_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator-=(const m_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator*=(const m_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator/=(const m_type s)
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

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator&=(const dtVector3 &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator&=(const dtVector<m_row, m_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> &dtVector3<m_type, m_row>::operator&=(const dtMatrix<m_row, 1, m_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtCommaInit<m_row, m_type> dtVector3<m_type, m_row>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator-() const
{
    return dtVector3(
        -m_elem[0],
        -m_elem[1],
        -m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator+(const dtVector3 &v) const
{
    return dtVector3(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator-(const dtVector3 &v) const
{
    return dtVector3(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator*(const dtVector3 &v) const
{
    return dtVector3(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator/(const dtVector3 &v) const
{
    m_type den[3];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];

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

    return dtVector3(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator+(const dtVector<m_row, m_type> &v) const
{
    return dtVector3(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator-(const dtVector<m_row, m_type> &v) const
{
    return dtVector3(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator*(const dtVector<m_row, m_type> &v) const
{
    return dtVector3(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator/(const dtVector<m_row, m_type> &v) const
{
    m_type den[3];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];

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

    return dtVector3(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator+(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector3(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator-(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector3(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator*(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector3(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator/(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type den[3];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];

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

    return dtVector3(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator+(const m_type s) const
{
    return dtVector3(
        m_elem[0] + s,
        m_elem[1] + s,
        m_elem[2] + s);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator-(const m_type s) const
{
    return dtVector3(
        m_elem[0] - s,
        m_elem[1] - s,
        m_elem[2] - s);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator*(const m_type s) const
{
    return dtVector3(
        m_elem[0] * s,
        m_elem[1] * s,
        m_elem[2] * s);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator/(const m_type s) const
{
    m_type den = s;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }

    return dtVector3(
        m_elem[0] / den,
        m_elem[1] / den,
        m_elem[2] / den);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator&(const dtVector3 &v) const
{
    return dtVector3<m_type, m_row>(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1],
        m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2],
        m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator&(const dtVector<m_row, m_type> &v) const
{
    return dtVector3<m_type, m_row>(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1],
        m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2],
        m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, m_row> dtVector3<m_type, m_row>::operator&(const dtMatrix<m_row, 1, m_type> &v) const
{
    return dtVector3<m_type, m_row>(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1],
        m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2],
        m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0]);
}

template <typename m_type, uint16_t m_row>
inline dtMatrix3<m_type, m_row, m_row> dtVector3<m_type, m_row>::operator&(const dtMatrix3<m_type, m_row, m_row> &m) const
{ // [v]x * Mat3, []x is skew-symmetric matrix
    return dtMatrix3<m_type, m_row, m_row>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <typename m_type, uint16_t m_row>
inline dtRotation<m_type, m_row, m_row> dtVector3<m_type, m_row>::operator&(const dtRotation<m_type, m_row, m_row> &m) const
{ // [v]x * RotMat, []x is skew-symmetric matrix
    return dtMatrix3<m_type, m_row, m_row>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <typename m_type, uint16_t m_row>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtVector3<m_type, m_row>::operator*(const dtMatrix<1, col, m_type> &m) const
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
inline m_type dtVector3<m_type, m_row>::dot(const dtVector3 &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector3<m_type, m_row>::dot(const dtVector<m_row, m_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtVector3<m_type, m_row>::dot(const dtMatrix<m_row, 1, m_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator==(const dtVector3 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator!=(const dtVector3 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator==(const dtVector<m_row, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator!=(const dtVector<m_row, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator==(const dtMatrix<m_row, 1, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtVector3<m_type, m_row>::operator!=(const dtMatrix<m_row, 1, m_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::Print(const char endChar)
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

//-- Private Member Function ------------------------------------------------//
template <typename m_type, uint16_t m_row>
inline void dtVector3<m_type, m_row>::CrossProduct(const m_type *v)
{
    m_type elem[m_row];

    elem[0] = m_elem[1] * v[2] - m_elem[2] * v[1];
    elem[1] = m_elem[2] * v[0] - m_elem[0] * v[2];
    elem[2] = m_elem[0] * v[1] - m_elem[1] * v[0];

    m_elem[0] = elem[0];
    m_elem[1] = elem[1];
    m_elem[2] = elem[2];
}

//-- Template Function ------------------------------------------------------//
// scalar * vector
template <typename type, uint16_t row>
inline dtVector3<type, row> operator+(const type s, const dtVector3<type, row> &v)
{
    return dtVector3<type, row>(
        v.m_elem[0] + s,
        v.m_elem[1] + s,
        v.m_elem[2] + s);
}

template <typename type, uint16_t row>
inline dtVector3<type, row> operator-(const type s, const dtVector3<type, row> &v)
{
    return dtVector3<type, row>(
        s - v.m_elem[0],
        s - v.m_elem[1],
        s - v.m_elem[2]);
}

template <typename type, uint16_t row>
inline dtVector3<type, row> operator*(const type s, const dtVector3<type, row> &v)
{
    return dtVector3<type, row>(
        v.m_elem[0] * s,
        v.m_elem[1] * s,
        v.m_elem[2] * s);
}

template <typename type, uint16_t row>
inline dtVector3<type, row> operator/(const type s, const dtVector3<type, row> &v)
{
    type den[3];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];

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

    return dtVector3<type, row>(
        s / den[0],
        s / den[1],
        s / den[2]);
}

typedef dtVector3<> dtVec3;

} // namespace dtMath

#endif // DTMATH_DTVECTOR3_TPP_