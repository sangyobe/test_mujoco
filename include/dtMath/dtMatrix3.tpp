/*!
\file       dtMatrix3.h
\brief      dtMath, 3x3 Matrix class, lighter and faster than general matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX3_TPP_
#define DTMATH_DTMATRIX3_TPP_

#include "dtMatrix3.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3()
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const m_type *element)
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const m_type *element, const size_t n_byte)
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const char c, const m_type *element, const size_t n_byte)
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
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const dtMatrix3 &m)
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const dtRotation<m_type, m_row, m_col> &m)
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
inline dtMatrix3<m_type, m_row, m_col>::dtMatrix3(const dtMatrix<m_row, m_col, m_type> &m)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetZero()
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
inline void dtMatrix3<m_type, m_row, m_col>::SetIdentity()
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
inline void dtMatrix3<m_type, m_row, m_col>::SetDiagonal(const m_type d1, const m_type d2, const m_type d3)
{
    m_elem[0] = d1;
    m_elem[4] = d2;
    m_elem[8] = d3;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetDiagonal(const m_type *element, const size_t n_byte)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetDiagonal(const dtVector<m_row, m_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetDiagonal(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetFill(const m_type value)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetElement(const m_type *element, const size_t n_byte)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetElement(
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
inline void dtMatrix3<m_type, m_row, m_col>::SetElement(const dtMatrix3 &m)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetElement(const dtRotation<m_type, m_row, m_col> &m)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetElement(const dtMatrix<m_row, m_col, m_type> &m)
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
template <uint16_t col>
inline void dtMatrix3<m_type, m_row, m_col>::SetRowVec(const uint16_t idxRow, const dtVector<col, m_type> &v)
{
    uint16_t maxCol = (m_col < col) ? m_col : col;

    for (uint16_t icol = 0; icol < maxCol; icol++)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetRowVec(const uint16_t idxRow, const dtVector3<m_type, m_col> &v)
{
    m_elem[idxRow * m_col] = v.m_elem[0];
    m_elem[idxRow * m_col + 1] = v.m_elem[1];
    m_elem[idxRow * m_col + 2] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte)
{
    uint16_t col = n_byte / sizeof(m_type);
    uint16_t maxCol = (m_col < col) ? m_col : col;

    for (uint16_t icol = 0; icol < maxCol; icol++)
    {
        m_elem[idxRow * m_col + icol] = v[icol];
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
template <uint16_t row>
inline void dtMatrix3<m_type, m_row, m_col>::SetColVec(const uint16_t idxCol, const dtVector<row, m_type> &v)
{
    uint16_t maxRow = (m_row < row) ? m_row : row;

    for (uint16_t irow = 0; irow < maxRow; irow++)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetColVec(const uint16_t idxCol, const dtVector3<m_type, m_row> &v)
{
    m_elem[idxCol] = v.m_elem[0];
    m_elem[m_col + idxCol] = v.m_elem[1];
    m_elem[m_col * 2 + idxCol] = v.m_elem[2];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte)
{
    uint16_t row = n_byte / sizeof(m_type);
    uint16_t maxRow = (m_row < row) ? m_row : row;

    for (uint16_t irow = 0; irow < maxRow; irow++)
    {
        m_elem[irow * m_col + idxCol] = v[irow];
    }
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtMatrix3<m_type, m_row, m_col>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
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
inline void dtMatrix3<m_type, m_row, m_col>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
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
inline const m_type *const dtMatrix3<m_type, m_row, m_col>::GetElementsAddr() const
{
    return m_elem;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, m_col> dtMatrix3<m_type, m_row, m_col>::GetRowVec(const uint16_t idxRow) const
{
    return dtVector3<m_type, m_row>(
        m_elem[idxRow * m_col],
        m_elem[idxRow * m_col + 1],
        m_elem[idxRow * m_col + 2]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, m_row> dtMatrix3<m_type, m_row, m_col>::GetColVec(const uint16_t idxCol) const
{
    return dtVector3<m_type, m_row>(
        m_elem[idxCol],
        m_elem[1 * m_col + idxCol],
        m_elem[2 * m_col + idxCol]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline int8_t dtMatrix3<m_type, m_row, m_col>::GetRowVec(const uint16_t idxRow, dtVector3<m_type, m_col> &v) const
{
    v.m_elem[0] = m_elem[idxRow * m_col];
    v.m_elem[1] = m_elem[idxRow * m_col + 1];
    v.m_elem[2] = m_elem[idxRow * m_col + 2];

    return 0;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline int8_t dtMatrix3<m_type, m_row, m_col>::GetColVec(const uint16_t idxCol, dtVector3<m_type, m_row> &v) const
{
    v.m_elem[0] = m_elem[idxCol];
    v.m_elem[1] = m_elem[1 * m_col + idxCol];
    v.m_elem[2] = m_elem[2 * m_col + idxCol];

    return 0;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::Transpose() const
{
    return dtMatrix3(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline m_type dtMatrix3<m_type, m_row, m_col>::Trace() const
{
    return (m_elem[0] + m_elem[4] + m_elem[8]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline m_type dtMatrix3<m_type, m_row, m_col>::GetNorm() const
{
    m_type sqSum =
        m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] + m_elem[4] * m_elem[4] + m_elem[5] * m_elem[5] +
        m_elem[6] * m_elem[6] + m_elem[7] * m_elem[7] + m_elem[8] * m_elem[8];

    return std::sqrt(sqSum);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline m_type dtMatrix3<m_type, m_row, m_col>::GetSqNorm() const
{
    m_type sqSum =
        m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] + m_elem[4] * m_elem[4] + m_elem[5] * m_elem[5] +
        m_elem[6] * m_elem[6] + m_elem[7] * m_elem[7] + m_elem[8] * m_elem[8];

    return sqSum;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtNoPivLU<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::NoPivLU() const
{
    return dtNoPivLU<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtPartialPivLU<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::PartialPivLU() const
{
    return dtPartialPivLU<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtLLT<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::LLT() const
{
    return dtLLT<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtLDLT<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::LDLT() const
{
    return dtLDLT<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtQR<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::QR() const
{
    return dtQR<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtSVD<m_row, m_col, m_type> dtMatrix3<m_type, m_row, m_col>::SVD() const
{
    return dtSVD<m_row, m_col, m_type>(*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::Inv(int8_t *isOk) const
{
    return dtMatrix3<m_type, m_row, m_col>(dtPartialPivLU<m_row, m_col, m_type>(*this).InverseArray(isOk));
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::PInv(int8_t *isOk, m_type tolerance) const
{
    return dtMatrix3<m_type, m_row, m_col>(dtSVD<m_row, m_col, m_type>(*this).InverseArray(isOk, tolerance));
}

/* Assignment operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator=(const dtMatrix3 &m)
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
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator+=(const dtMatrix3 &m)
{
    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator-=(const dtMatrix3 &m)
{
    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator=(const dtRotation<m_type, m_row, m_col> &m)
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
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator+=(const dtRotation<m_type, m_row, m_col> &m)
{
    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator-=(const dtRotation<m_type, m_row, m_col> &m)
{
    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator=(const dtMatrix<m_row, m_col, m_type> &m)
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
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator+=(const dtMatrix<m_row, m_col, m_type> &m)
{
    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator-=(const dtMatrix<m_row, m_col, m_type> &m)
{
    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator=(const m_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;
    m_elem[3] = s;
    m_elem[4] = s;
    m_elem[5] = s;
    m_elem[6] = s;
    m_elem[7] = s;
    m_elem[8] = s;

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator+=(const m_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;
    m_elem[3] += s;
    m_elem[4] += s;
    m_elem[5] += s;
    m_elem[6] += s;
    m_elem[7] += s;
    m_elem[8] += s;

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator-=(const m_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;
    m_elem[3] -= s;
    m_elem[4] -= s;
    m_elem[5] -= s;
    m_elem[6] -= s;
    m_elem[7] -= s;
    m_elem[8] -= s;

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator*=(const m_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;
    m_elem[4] *= s;
    m_elem[5] *= s;
    m_elem[6] *= s;
    m_elem[7] *= s;
    m_elem[8] *= s;

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> &dtMatrix3<m_type, m_row, m_col>::operator/=(const m_type s)
{
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    m_elem[0] /= scalar;
    m_elem[1] /= scalar;
    m_elem[2] /= scalar;
    m_elem[3] /= scalar;
    m_elem[4] /= scalar;
    m_elem[5] /= scalar;
    m_elem[6] /= scalar;
    m_elem[7] /= scalar;
    m_elem[8] /= scalar;

    return (*this);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtCommaInit<m_row * m_col, m_type> dtMatrix3<m_type, m_row, m_col>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row * m_col, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator-() const
{
    return dtMatrix3(*this) *= -1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator+(const dtMatrix3 &m) const
{
    return dtMatrix3(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator-(const dtMatrix3 &m) const
{
    return dtMatrix3(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator+(const dtRotation<m_type, m_row, m_col> &m) const
{
    return dtMatrix3(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator-(const dtRotation<m_type, m_row, m_col> &m) const
{
    return dtMatrix3(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator+(const dtMatrix<m_row, m_col, m_type> &m) const
{
    return dtMatrix3(*this) += m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator-(const dtMatrix<m_row, m_col, m_type> &m) const
{
    return dtMatrix3(*this) -= m;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator+(const m_type s) const
{
    return dtMatrix3(*this) += s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator-(const m_type s) const
{
    return dtMatrix3(*this) -= s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator*(const m_type s) const
{
    return dtMatrix3(*this) *= s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator/(const m_type s) const
{
    return dtMatrix3(*this) /= s;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtMatrix3<m_type, m_row, m_col>::operator*(const dtMatrix<m_col, col, m_type> &m) const
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
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator*(const dtMatrix3 &m) const
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
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator*(const dtRotation<m_type, m_row, m_col> &m) const
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
inline dtVector<m_row, m_type> dtMatrix3<m_type, m_row, m_col>::operator*(const dtVector<m_col, m_type> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return dtVector<m_row, m_type>(vec);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, m_row> dtMatrix3<m_type, m_row, m_col>::operator*(const dtVector3<m_type, m_col> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return dtVector3<m_type, m_row>(vec);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator&(const dtVector<m_col, m_type> &v) const
{ // Mat3 * [v]x, []x is skew-symmetric matrix
    return dtMatrix3(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix3<m_type, m_row, m_col> dtMatrix3<m_type, m_row, m_col>::operator&(const dtVector3<m_type, m_col> &v) const
{ // Mat3 * [v]x, []x is skew-symmetric matrix
    return dtMatrix3(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtMatrix3<m_type, m_row, m_col>::operator==(const dtMatrix3 &m) const
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
inline bool dtMatrix3<m_type, m_row, m_col>::operator!=(const dtMatrix3 &m) const
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
inline bool dtMatrix3<m_type, m_row, m_col>::operator==(const dtRotation<m_type, m_row, m_col> &m) const
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
inline bool dtMatrix3<m_type, m_row, m_col>::operator!=(const dtRotation<m_type, m_row, m_col> &m) const
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
inline bool dtMatrix3<m_type, m_row, m_col>::operator==(const dtMatrix<m_row, m_col, m_type> &m) const
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
inline bool dtMatrix3<m_type, m_row, m_col>::operator!=(const dtMatrix<m_row, m_col, m_type> &m) const
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
inline void dtMatrix3<m_type, m_row, m_col>::Print(const char endChar)
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

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <typename type, uint16_t row>
inline dtMatrix3<type, row> operator*(const type s, const dtMatrix3<type, row> &m)
{
    return dtMatrix3<type, row>(m) *= s;
}

typedef dtMatrix3<double> dtMat3d;
typedef dtMatrix3<float> dtMat3;

} // namespace dtMath

#endif // DTMATH_DTMATRIX3_TPP_
