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

#ifndef DTMATH_DTVECTOR_TPP_
#define DTMATH_DTVECTOR_TPP_

#include "dtVector.h"

namespace dtMath
{

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type>::dtVector()
{
    memset(m_elem, 0, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type>::dtVector(const m_type *element)
{
    memcpy(m_elem, element, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type>::dtVector(const m_type *element, const size_t n_byte)
{
    size_t vecSz = sizeof(m_type) * m_row;

    if (vecSz > n_byte)
    {
        memset(m_elem, 0, sizeof(m_type) * m_row);
        memcpy(m_elem, element, n_byte);
    }
    else
        memcpy(m_elem, element, vecSz);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type>::dtVector(const dtVector &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type>::dtVector(const dtMatrix<m_row, 1, m_type> &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetZero()
{
    memset(m_elem, 0, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetFill(const m_type value)
{
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = value;
        m_elem[irow + 1] = value;
        m_elem[irow + 2] = value;
        m_elem[irow + 3] = value;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] = value;
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetElement(const m_type *element, const size_t n_byte)
{
    size_t vecSz = sizeof(m_type) * m_row;

    if (vecSz > n_byte)
        memcpy(m_elem, element, n_byte);
    else
        memcpy(m_elem, element, vecSz);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetElement(const dtVector &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetElement(const dtMatrix<m_row, 1, m_type> &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);
}

// template<uint16_t m_row, typename m_type>
// inline void dtVector<m_row, m_type>::SetElement(const m_type elem0, ...)
//{
//     va_list ap;
//     va_start(ap, elem0);
//
//     m_elem[0] = elem0;
//
//     for (uint16_t irow = 1; irow < m_row; ++irow)
//     {
//         m_elem[irow] = (m_type)(va_arg(ap, double));
//     }
//
//     va_end(ap);
// }

template <uint16_t m_row, typename m_type>
template <uint16_t row>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const dtVector<row, m_type> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row)
        rowSz = row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(m_type) * rowSz);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const m_type *v, const size_t n_byte)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    uint16_t row = (uint16_t)(n_byte / sizeof(m_type));
    if (rowSz > row)
        rowSz = row;

    memcpy(&m_elem[idxRow], v, sizeof(m_type) * rowSz);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const dtVector3<m_type, 3> &v)
{
    if (idxRow >= m_row)
        return;

    switch (m_row - idxRow)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    }
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const dtVector4<m_type, 4> &v)
{
    if (idxRow >= m_row)
        return;

    switch (m_row - idxRow)
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
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    }
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const dtVector6<m_type, 6> &v)
{
    if (idxRow >= m_row)
        return;

    switch (m_row - idxRow)
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
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        m_elem[idxRow + 4] = v.m_elem[4];
        m_elem[idxRow + 5] = v.m_elem[5];
        break;
    }
}

template <uint16_t m_row, typename m_type>
template <uint16_t row>
inline void dtVector<m_row, m_type>::SetBlock(const uint16_t idxRow, const dtMatrix<row, 1, m_type> &v)
{
    if (idxRow >= m_row)
        return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row)
        rowSz = row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(m_type) * rowSz);
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetSwap(const uint16_t i, const uint16_t j)
{
    m_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::SetNormalize()
{
    uint16_t cnt = 0;
    uint16_t irow = 0;
    m_type norm = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    for (cnt = m_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= norm;
        m_elem[irow + 1] /= norm;
        m_elem[irow + 2] /= norm;
        m_elem[irow + 3] /= norm;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] /= norm;
}

template <uint16_t m_row, typename m_type>
inline const m_type *const dtVector<m_row, m_type>::GetElementsAddr() const
{
    return m_elem;
}

template <uint16_t m_row, typename m_type>
template <uint16_t row>
inline dtVector<row, m_type> dtVector<m_row, m_type>::GetBlock(const uint16_t idx)
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

template <uint16_t m_row, typename m_type>
inline dtVector3<m_type, 3> dtVector<m_row, m_type>::GetBlockVec3(const uint16_t idx)
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

template <uint16_t m_row, typename m_type>
inline dtVector4<m_type, 4> dtVector<m_row, m_type>::GetBlockVec4(const uint16_t idx)
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

template <uint16_t m_row, typename m_type>
inline dtVector6<m_type, 6> dtVector<m_row, m_type>::GetBlockVec6(const uint16_t idx)
{
    m_type elem[6] = {
        0,
    };

    if (idx >= m_row)
        return dtVector6<m_type, 6>(elem);

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
    case 4:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        break;
    case 5:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        elem[4] = m_elem[idx + 4];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        elem[4] = m_elem[idx + 4];
        elem[5] = m_elem[idx + 5];
    };

    return dtVector6<m_type, 6>(elem);
}

template <uint16_t m_row, typename m_type>
template <uint16_t row>
inline int8_t dtVector<m_row, m_type>::GetBlock(const uint16_t idx, dtVector<row, m_type> &v)
{
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row)
        return -1;
    if (rowSize > row)
        rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(m_type) * rowSize);

    return 0;
}

template <uint16_t m_row, typename m_type>
inline int8_t dtVector<m_row, m_type>::GetBlockVec3(const uint16_t idx, dtVector3<m_type, 3> &v)
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

template <uint16_t m_row, typename m_type>
inline int8_t dtVector<m_row, m_type>::GetBlockVec4(const uint16_t idx, dtVector4<m_type, 4> &v)
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

template <uint16_t m_row, typename m_type>
inline int8_t dtVector<m_row, m_type>::GetBlockVec6(const uint16_t idx, dtVector6<m_type, 6> &v)
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
    case 4:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        break;
    case 5:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        v.m_elem[4] = m_elem[idx + 4];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        v.m_elem[4] = m_elem[idx + 4];
        v.m_elem[5] = m_elem[idx + 5];
    };

    return 0;
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::GetNorm() const
{
    m_type rtn = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        rtn += m_elem[irow] * m_elem[irow];
        rtn += m_elem[irow + 1] * m_elem[irow + 1];
        rtn += m_elem[irow + 2] * m_elem[irow + 2];
        rtn += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        rtn += m_elem[irow] * m_elem[irow];
    }

    return std::sqrt(rtn);
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::GetSqNorm() const
{
    m_type rtn = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        rtn += m_elem[irow] * m_elem[irow];
        rtn += m_elem[irow + 1] * m_elem[irow + 1];
        rtn += m_elem[irow + 2] * m_elem[irow + 2];
        rtn += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        rtn += m_elem[irow] * m_elem[irow];
    }

    return rtn;
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::GetSum() const
{
    m_type rtn = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        rtn += m_elem[irow];
        rtn += m_elem[irow + 1];
        rtn += m_elem[irow + 2];
        rtn += m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        rtn += m_elem[irow];
    }

    return rtn;
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::GetNormalized() const
{
    m_type vec[m_row];
    uint16_t cnt = 0;
    uint16_t irow = 0;
    m_type norm = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    for (cnt = m_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] / norm;
        vec[irow + 1] = m_elem[irow + 1] / norm;
        vec[irow + 2] = m_elem[irow + 2] / norm;
        vec[irow + 3] = m_elem[irow + 3] / norm;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        vec[irow] = m_elem[irow] / norm;

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, typename m_type>
inline dtMatrix<3, 3, m_type> dtVector<m_row, m_type>::GetSkew() const
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");

    return dtMatrix<3, 3, m_type>(
        0, -m_elem[2], m_elem[1],
        m_elem[2], 0, -m_elem[0],
        -m_elem[1], m_elem[0], 0);
}

template <uint16_t m_row, typename m_type>
inline dtMatrix<1, m_row, m_type> dtVector<m_row, m_type>::Transpose() const
{
    return dtMatrix<1, m_row, m_type>(m_elem);
}

/* Assignment operators */
template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const dtVector &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const dtVector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const dtVector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const dtVector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const dtVector &v)
{
    uint16_t cnt, irow = 0;
    m_type den[m_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 1] < 0)
                den[irow + 1] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 1] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 2] < 0)
                den[irow + 2] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 2] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 3] < 0)
                den[irow + 3] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 3] = std::numeric_limits<m_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
        m_elem[irow + 1] /= den[irow + 1];
        m_elem[irow + 2] /= den[irow + 2];
        m_elem[irow + 3] /= den[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const dtVector3<m_type, m_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const dtVector3<m_type, m_row> &v)
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

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const dtVector4<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const dtVector4<m_type, m_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const dtVector4<m_type, m_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const dtVector4<m_type, m_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const dtVector4<m_type, m_row> &v)
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

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const dtVector6<m_type, m_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const dtVector6<m_type, m_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const dtVector6<m_type, m_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const dtVector6<m_type, m_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const dtVector6<m_type, m_row> &v)
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

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const dtMatrix<m_row, 1, m_type> &v)
{
    memcpy(m_elem, v.m_elem, sizeof(m_type) * m_row);

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const dtMatrix<m_row, 1, m_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u, irow; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const dtMatrix<m_row, 1, m_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const dtMatrix<m_row, 1, m_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u, irow; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const dtMatrix<m_row, 1, m_type> &v)
{
    uint16_t cnt, irow = 0;
    m_type den[m_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = m_row >> 2u, irow; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 1] < 0)
                den[irow + 1] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 1] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 2] < 0)
                den[irow + 2] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 2] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 3] < 0)
                den[irow + 3] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 3] = std::numeric_limits<m_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
        m_elem[irow + 1] /= den[irow + 1];
        m_elem[irow + 2] /= den[irow + 2];
        m_elem[irow + 3] /= den[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator=(const m_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = s;
        m_elem[irow + 1] = s;
        m_elem[irow + 2] = s;
        m_elem[irow + 3] = s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] = s;
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator+=(const m_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += s;
        m_elem[irow + 1] += s;
        m_elem[irow + 2] += s;
        m_elem[irow + 3] += s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += s;
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator-=(const m_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= s;
        m_elem[irow + 1] /= s;
        m_elem[irow + 2] /= s;
        m_elem[irow + 3] /= s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] /= s;
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator*=(const m_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= s;
        m_elem[irow + 1] *= s;
        m_elem[irow + 2] *= s;
        m_elem[irow + 3] *= s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= s;
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator/=(const m_type s)
{
    m_type scalar = s;
    uint16_t cnt, irow = 0;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= scalar;
        m_elem[irow + 1] /= scalar;
        m_elem[irow + 2] /= scalar;
        m_elem[irow + 3] /= scalar;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] /= scalar;
    }

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator&=(const dtVector3<m_type, m_row> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator&=(const dtVector<m_row, m_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> &dtVector<m_row, m_type>::operator&=(const dtMatrix<m_row, 1, m_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t m_row, typename m_type>
inline dtCommaInit<m_row, m_type> dtVector<m_row, m_type>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row, m_type>(m_elem);
}

/* Arithmetic operators */
template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-() const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = -m_elem[irow];
        vec[irow + 1] = -m_elem[irow + 1];
        vec[irow + 2] = -m_elem[irow + 2];
        vec[irow + 3] = -m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = -m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const dtVector &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const dtVector &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const dtVector &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const dtVector &v) const
{
    m_type vec[m_row];
    m_type den[m_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 1] < 0)
                den[irow + 1] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 1] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 2] < 0)
                den[irow + 2] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 2] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 3] < 0)
                den[irow + 3] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 3] = std::numeric_limits<m_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
        vec[irow + 1] = m_elem[irow + 1] / den[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] / den[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] / den[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const dtVector3<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const dtVector3<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const dtVector3<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const dtVector3<m_type, m_row> &v) const
{
    m_type vec[m_row];
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const dtVector4<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];
    vec[3] = m_elem[3] + v.m_elem[3];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const dtVector4<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];
    vec[3] = m_elem[3] - v.m_elem[3];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const dtVector4<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];
    vec[3] = m_elem[3] * v.m_elem[3];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const dtVector4<m_type, m_row> &v) const
{
    m_type vec[m_row];
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[3] = m_elem[3] / den;

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const dtVector6<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];
    vec[3] = m_elem[3] + v.m_elem[3];
    vec[4] = m_elem[4] + v.m_elem[4];
    vec[5] = m_elem[5] + v.m_elem[5];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const dtVector6<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];
    vec[3] = m_elem[3] - v.m_elem[3];
    vec[4] = m_elem[4] - v.m_elem[4];
    vec[5] = m_elem[5] - v.m_elem[5];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const dtVector6<m_type, m_row> &v) const
{
    m_type vec[m_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];
    vec[3] = m_elem[3] * v.m_elem[3];
    vec[4] = m_elem[4] * v.m_elem[4];
    vec[5] = m_elem[5] * v.m_elem[5];

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const dtVector6<m_type, m_row> &v) const
{
    m_type vec[m_row];
    m_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[3] = m_elem[3] / den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[4] = m_elem[4] / den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }
    vec[5] = m_elem[5] / den;

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type vec[m_row];
    m_type den[m_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 1] < 0)
                den[irow + 1] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 1] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 2] < 0)
                den[irow + 2] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 2] = std::numeric_limits<m_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow + 3] < 0)
                den[irow + 3] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow + 3] = std::numeric_limits<m_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
        vec[irow + 1] = m_elem[irow + 1] / den[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] / den[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] / den[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<m_type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<m_type>::epsilon();
            else
                den[irow] = std::numeric_limits<m_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator+(const m_type s) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + s;
        vec[irow + 1] = m_elem[irow + 1] + s;
        vec[irow + 2] = m_elem[irow + 2] + s;
        vec[irow + 3] = m_elem[irow + 3] + s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + s;
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator-(const m_type s) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - s;
        vec[irow + 1] = m_elem[irow + 1] - s;
        vec[irow + 2] = m_elem[irow + 2] - s;
        vec[irow + 3] = m_elem[irow + 3] - s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - s;
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator*(const m_type s) const
{
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * s;
        vec[irow + 1] = m_elem[irow + 1] * s;
        vec[irow + 2] = m_elem[irow + 2] * s;
        vec[irow + 3] = m_elem[irow + 3] * s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * s;
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator/(const m_type s) const
{
    m_type den = s;
    m_type vec[m_row];
    uint16_t cnt, irow = 0;

    if (std::abs(den) < std::numeric_limits<m_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<m_type>::epsilon();
        else
            den = std::numeric_limits<m_type>::epsilon();
    }

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] / den;
        vec[irow + 1] = m_elem[irow + 1] / den;
        vec[irow + 2] = m_elem[irow + 2] / den;
        vec[irow + 3] = m_elem[irow + 3] / den;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] / den;
    }

    return dtVector(vec);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator&(const dtVector3<m_type, m_row> &v) const
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");

    m_type elem[m_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return dtVector<m_row, m_type>(elem);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator&(const dtVector<m_row, m_type> &v) const
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");

    m_type elem[m_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return dtVector<m_row, m_type>(elem);
}

template <uint16_t m_row, typename m_type>
inline dtVector<m_row, m_type> dtVector<m_row, m_type>::operator&(const dtMatrix<m_row, 1, m_type> &v) const
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");

    m_type elem[m_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return dtVector<m_row, m_type>(elem);
}

template <uint16_t m_row, typename m_type>
inline dtMatrix3<m_type, m_row, m_row> dtVector<m_row, m_type>::operator&(const dtMatrix3<m_type, m_row, m_row> &m) const
{ // [v]x * Mat3, []x is skew-symmetric matrix
    return dtMatrix3<m_type, m_row, m_row>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <uint16_t m_row, typename m_type>
inline dtRotation<m_type, m_row, m_row> dtVector<m_row, m_type>::operator&(const dtRotation<m_type, m_row, m_row> &m) const
{ // [v]x * RotMat, []x is skew-symmetric matrix
    return dtMatrix3<m_type, m_row, m_row>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <uint16_t m_row, typename m_type>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtVector<m_row, m_type>::operator*(const dtMatrix<1, col, m_type> &m) const
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

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::dot(const dtVector &v) const
{
    m_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::dot(const dtVector3<m_type, m_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::dot(const dtVector4<m_type, m_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::dot(const dtVector6<m_type, m_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

template <uint16_t m_row, typename m_type>
inline m_type dtVector<m_row, m_type>::dot(const dtMatrix<m_row, 1, m_type> &v) const
{
    m_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

/* Comparison operators */
template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator==(const dtVector &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance)
            return false;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return false;
    }

    return true;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator!=(const dtVector &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance)
            return true;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return true;
    }

    return false;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator==(const dtVector3<m_type, m_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;

    return true;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator!=(const dtVector3<m_type, m_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;

    return false;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator==(const dtVector4<m_type, m_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return false;

    return true;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator!=(const dtVector4<m_type, m_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance)
        return true;

    return false;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator==(const dtVector6<m_type, m_row> &v) const
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

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator!=(const dtVector6<m_type, m_row> &v) const
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

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator==(const dtMatrix<m_row, 1, m_type> &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance)
            return false;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return false;
    }

    return true;
}

template <uint16_t m_row, typename m_type>
inline bool dtVector<m_row, m_type>::operator!=(const dtMatrix<m_row, 1, m_type> &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance)
            return true;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance)
            return true;
    }

    return false;
}

template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::Print(const char endChar)
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
template <uint16_t m_row, typename m_type>
inline void dtVector<m_row, m_type>::CrossProduct(const m_type *v)
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");

    m_type elem[m_row];

    elem[0] = m_elem[1] * v[2] - m_elem[2] * v[1];
    elem[1] = m_elem[2] * v[0] - m_elem[0] * v[2];
    elem[2] = m_elem[0] * v[1] - m_elem[1] * v[0];

    m_elem[0] = elem[0];
    m_elem[1] = elem[1];
    m_elem[2] = elem[2];
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <uint16_t row, typename type>
inline dtVector<row, type> operator+(const type s, const dtVector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = v.m_elem[irow] + s;
        vec[irow + 1] = v.m_elem[irow + 1] + s;
        vec[irow + 2] = v.m_elem[irow + 2] + s;
        vec[irow + 3] = v.m_elem[irow + 3] + s;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = v.m_elem[irow] + s;
    }

    return dtVector<row, type>(vec);
}

// scalar - vector
template <uint16_t row, typename type>
inline dtVector<row, type> operator-(const type s, const dtVector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = s - v.m_elem[irow];
        vec[irow + 1] = s - v.m_elem[irow + 1];
        vec[irow + 2] = s - v.m_elem[irow + 2];
        vec[irow + 3] = s - v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = s - v.m_elem[irow];
    }

    return dtVector<row, type>(vec);
}

// scalar * vector
template <uint16_t row, typename type>
inline dtVector<row, type> operator*(const type s, const dtVector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = v.m_elem[irow] * s;
        vec[irow + 1] = v.m_elem[irow + 1] * s;
        vec[irow + 2] = v.m_elem[irow + 2] * s;
        vec[irow + 3] = v.m_elem[irow + 3] * s;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = v.m_elem[irow] * s;
    }

    return dtVector<row, type>(vec);
}

// scalar / vector
template <uint16_t row, typename type>
inline dtVector<row, type> operator/(const type s, const dtVector<row, type> &v)
{
    type vec[row];
    type den[row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<type>::epsilon();
            else
                den[irow] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 1] < 0)
                den[irow + 1] = -std::numeric_limits<type>::epsilon();
            else
                den[irow + 1] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 2] < 0)
                den[irow + 2] = -std::numeric_limits<type>::epsilon();
            else
                den[irow + 2] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 3] < 0)
                den[irow + 3] = -std::numeric_limits<type>::epsilon();
            else
                den[irow + 3] = std::numeric_limits<type>::epsilon();
        }
        vec[irow] = s / den[irow];
        vec[irow + 1] = s / den[irow + 1];
        vec[irow + 2] = s / den[irow + 2];
        vec[irow + 3] = s / den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow] < 0)
                den[irow] = -std::numeric_limits<type>::epsilon();
            else
                den[irow] = std::numeric_limits<type>::epsilon();
        }
        vec[irow] = s / den[irow];
    }

    return dtVector<row, type>(vec);
}

} // namespace dtMath

#endif // DTMATH_DTVECTOR_TPP_