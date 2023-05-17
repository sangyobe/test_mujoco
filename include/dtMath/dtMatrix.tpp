/*!
\file       dtMatrix.h
\brief      dtMath, General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX_TPP_
#define DTMATH_DTMATRIX_TPP_

#include "dtMatrix.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type>::dtMatrix()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type>::dtMatrix(const m_type *element)
{
    memcpy(m_elem, element, sizeof(m_type) * m_row * m_col);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type>::dtMatrix(const m_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(m_type) * m_row * m_col;

    if (matSz <= n_byte)
        memcpy(m_elem, element, matSz);
    else
    {
        memset(m_elem, 0, matSz);
        memcpy(m_elem, element, n_byte);
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type>::dtMatrix(const char c, const m_type *element, const size_t n_byte)
{
    if (c == 'a')
    {
        size_t matSz = sizeof(m_type) * m_row * m_col;

        if (matSz <= n_byte)
        {
            // memcpy(m_elem, element, matSz);
            uint16_t i, j, k;
            k = 0;
            for (i = 0; i < m_row; ++i)
                for (j = 0; j < m_col; ++j)
                    m_elem[i * m_col + j] = element[k++];
        }
        else
        {
            memset(m_elem, 0, matSz);
            memcpy(m_elem, element, n_byte);
        }
    }

    else if (c == 'd')
    {
        memset(m_elem, 0, sizeof(m_type) * m_row * m_col);

        uint16_t num = (m_row > m_col) ? m_col : m_row;
        uint16_t offset = m_col + 1;
        uint16_t elemNum = (uint16_t)(n_byte / sizeof(m_type));

        num = (num > elemNum) ? elemNum : num;

        for (uint16_t i = 0; i < num; i++)
            m_elem[i * offset] = element[i];
    }

    else
    {
        memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type>::dtMatrix(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetZero()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetIdentity()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;

    for (uint16_t i = 0; i < num; i++)
        m_elem[i * offset] = 1;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetDiagonal(const m_type *element, const size_t n_byte)
{
    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;
    uint16_t elemNum = (uint16_t)(n_byte / sizeof(m_type));

    num = (num > elemNum) ? elemNum : num;

    for (uint16_t i = 0; i < num; i++)
        m_elem[i * offset] = element[i];
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetFill(const m_type value)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] = value;
        m_elem[i + 1] = value;
        m_elem[i + 2] = value;
        m_elem[i + 3] = value;
    }

    for (cnt = (m_row * m_col) % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = value;
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetElement(const m_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(m_type) * m_row * m_col;

    if (matSz <= n_byte)
        memcpy(m_elem, element, matSz);
    else
        memcpy(m_elem, element, n_byte);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t row, uint16_t col>
inline void dtMatrix<m_row, m_col, m_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtMatrix<row, col, m_type> &m)
{
    for (uint16_t irow = 0; irow < row; ++irow)
        for (uint16_t icol = 0; icol < col; ++icol)
            m_elem[(irow + idxRow) * m_col + idxCol + icol] = m.m_elem[irow * col + icol];
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtMatrix3<m_type, 3, 3> &m)
{
    m_elem[idxRow * m_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * m_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * m_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * m_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * m_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * m_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * m_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * m_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * m_col + idxCol + 2] = m.m_elem[8];
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const dtRotation<m_type, 3, 3> &m)
{
    m_elem[idxRow * m_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * m_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * m_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * m_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * m_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * m_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * m_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * m_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * m_col + idxCol + 2] = m.m_elem[8];
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline void dtMatrix<m_row, m_col, m_type>::SetRowVec(const uint16_t idxRow, const dtVector<col, m_type> &v)
{
    uint16_t maxCol = (m_col < col) ? m_col : col;
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = maxCol >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
        m_elem[idxRow * m_col + icol + 1] = v.m_elem[icol + 1];
        m_elem[idxRow * m_col + icol + 2] = v.m_elem[icol + 2];
        m_elem[idxRow * m_col + icol + 3] = v.m_elem[icol + 3];
    }

    for (cnt = maxCol % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetRowVec(const uint16_t idxRow, const dtVector3<m_type, 3> &v)
{
    if (m_col >= 3)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetRowVec(const uint16_t idxRow, const dtVector4<m_type, 4> &v)
{
    if (m_col >= 4)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
        m_elem[m_col * idxRow + 3] = v.m_elem[3];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetRowVec(const uint16_t idxRow, const dtVector6<m_type, 6> &v)
{
    if (m_col >= 6)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
        m_elem[m_col * idxRow + 3] = v.m_elem[3];
        m_elem[m_col * idxRow + 4] = v.m_elem[4];
        m_elem[m_col * idxRow + 5] = v.m_elem[5];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte)
{
    uint16_t col = (uint16_t)(n_byte / sizeof(m_type));
    uint16_t maxCol = (m_col < col) ? m_col : col;
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = maxCol >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = v[icol];
        m_elem[idxRow * m_col + icol + 1] = v[icol + 1];
        m_elem[idxRow * m_col + icol + 2] = v[icol + 2];
        m_elem[idxRow * m_col + icol + 3] = v[icol + 3];
    }

    for (cnt = maxCol % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = v[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t row>
inline void dtMatrix<m_row, m_col, m_type>::SetColVec(const uint16_t idxCol, const dtVector<row, m_type> &v)
{
    uint16_t maxRow = (m_row < row) ? m_row : row;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = maxRow >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
        m_elem[(irow + 1) * m_col + idxCol] = v.m_elem[irow + 1];
        m_elem[(irow + 2) * m_col + idxCol] = v.m_elem[irow + 2];
        m_elem[(irow + 3) * m_col + idxCol] = v.m_elem[irow + 3];
    }

    for (cnt = maxRow % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetColVec(const uint16_t idxCol, const dtVector3<m_type, 3> &v)
{
    if (m_row >= 3)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetColVec(const uint16_t idxCol, const dtVector4<m_type, 4> &v)
{
    if (m_row >= 4)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
        m_elem[m_col * 3 + idxCol] = v.m_elem[3];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetColVec(const uint16_t idxCol, const dtVector6<m_type, 6> &v)
{
    if (m_row >= 6)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
        m_elem[m_col * 3 + idxCol] = v.m_elem[3];
        m_elem[m_col * 4 + idxCol] = v.m_elem[4];
        m_elem[m_col * 5 + idxCol] = v.m_elem[5];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte)
{
    uint16_t row = n_byte / sizeof(m_type);
    uint16_t maxRow = (m_row < row) ? m_row : row;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = maxRow >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = v[irow];
        m_elem[(irow + 1) * m_col + idxCol] = v[irow + 1];
        m_elem[(irow + 2) * m_col + idxCol] = v[irow + 2];
        m_elem[(irow + 3) * m_col + idxCol] = v[irow + 3];
    }

    for (cnt = maxRow % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = v[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
{
    m_type tmpVec[m_col];
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = m_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        tmpVec[icol] = m_elem[idxRow1 * m_col + icol];
        m_elem[idxRow1 * m_col + icol] = m_elem[idxRow2 * m_col + icol];
        m_elem[idxRow2 * m_col + icol] = tmpVec[icol];

        tmpVec[icol + 1] = m_elem[idxRow1 * m_col + icol + 1];
        m_elem[idxRow1 * m_col + icol + 1] = m_elem[idxRow2 * m_col + icol + 1];
        m_elem[idxRow2 * m_col + icol + 1] = tmpVec[icol + 1];

        tmpVec[icol + 2] = m_elem[idxRow1 * m_col + icol + 2];
        m_elem[idxRow1 * m_col + icol + 2] = m_elem[idxRow2 * m_col + icol + 2];
        m_elem[idxRow2 * m_col + icol + 2] = tmpVec[icol + 2];

        tmpVec[icol + 3] = m_elem[idxRow1 * m_col + icol + 3];
        m_elem[idxRow1 * m_col + icol + 3] = m_elem[idxRow2 * m_col + icol + 3];
        m_elem[idxRow2 * m_col + icol + 3] = tmpVec[icol + 3];
    }

    for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
    {
        tmpVec[icol] = m_elem[idxRow1 * m_col + icol];
        m_elem[idxRow1 * m_col + icol] = m_elem[idxRow2 * m_col + icol];
        m_elem[idxRow2 * m_col + icol] = tmpVec[icol];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
{
    m_type tmpVec[m_row];
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        tmpVec[irow] = m_elem[irow * m_col + idxCol1];
        m_elem[irow * m_col + idxCol1] = m_elem[irow * m_col + idxCol2];
        m_elem[irow * m_col + idxCol2] = tmpVec[irow];

        tmpVec[irow + 1] = m_elem[(irow + 1) * m_col + idxCol1];
        m_elem[(irow + 1) * m_col + idxCol1] = m_elem[(irow + 1) * m_col + idxCol2];
        m_elem[(irow + 1) * m_col + idxCol2] = tmpVec[irow + 1];

        tmpVec[irow + 2] = m_elem[(irow + 2) * m_col + idxCol1];
        m_elem[(irow + 2) * m_col + idxCol1] = m_elem[(irow + 2) * m_col + idxCol2];
        m_elem[(irow + 2) * m_col + idxCol2] = tmpVec[irow + 2];

        tmpVec[irow + 3] = m_elem[(irow + 3) * m_col + idxCol1];
        m_elem[(irow + 3) * m_col + idxCol1] = m_elem[(irow + 3) * m_col + idxCol2];
        m_elem[(irow + 3) * m_col + idxCol2] = tmpVec[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        tmpVec[irow] = m_elem[irow * m_col + idxCol1];
        m_elem[irow * m_col + idxCol1] = m_elem[irow * m_col + idxCol2];
        m_elem[irow * m_col + idxCol2] = tmpVec[irow];
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline const m_type *const dtMatrix<m_row, m_col, m_type>::GetElementsAddr() const
{
    return m_elem;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t row, uint16_t col>
inline dtMatrix<row, col, m_type> dtMatrix<m_row, m_col, m_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol)
{
    uint16_t irow, icol, cnt;
    m_type elem[row * col] = {
        0,
    };
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;

    if (idxRow >= m_row)
        return dtMatrix<row, col, m_type>(elem);
    if (idxCol >= m_col)
        return dtMatrix<row, col, m_type>(elem);
    if (rowSize > row)
        rowSize = row;
    if (colSize > col)
        colSize = col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            elem[irow * col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            elem[irow * col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            elem[irow * col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return dtMatrix<row, col, m_type>(elem);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtMatrix<m_row, m_col, m_type>::GetRowVec(const uint16_t idxRow) const
{
    m_type vec[m_col];
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = m_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        vec[icol] = m_elem[idxRow * m_col + icol];
        vec[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        vec[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        vec[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
    {
        vec[icol] = m_elem[idxRow * m_col + icol];
    }

    return dtVector<m_col, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtMatrix<m_row, m_col, m_type>::GetColVec(const uint16_t idxCol) const
{
    m_type vec[m_row];
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow * m_col + idxCol];
        vec[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        vec[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        vec[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow * m_col + idxCol];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t row, uint16_t col>
inline int8_t dtMatrix<m_row, m_col, m_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, dtMatrix<row, col, m_type> &m)
{
    uint16_t irow, icol, cnt;
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;

    if (idxRow >= m_row)
        return -1;
    if (idxCol >= m_col)
        return -1;
    if (rowSize > row)
        rowSize = row;
    if (colSize > col)
        colSize = col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            m.m_elem[irow * col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            m.m_elem[irow * col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            m.m_elem[irow * col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtMatrix<m_row, m_col, m_type>::GetRowVec(const uint16_t idxRow, dtVector<m_col, m_type> &v) const
{
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = m_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
        v.m_elem[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        v.m_elem[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        v.m_elem[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtMatrix<m_row, m_col, m_type>::GetColVec(const uint16_t idxCol, dtVector<m_row, m_type> &v) const
{
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::GetCscMat() const
{
    return dtCscMatrix<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_col, m_row, m_type> dtMatrix<m_row, m_col, m_type>::Transpose() const
{
    m_type mat[m_row * m_col];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < m_row; ++irow)
    {
        for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            mat[irow + m_row * icol] = m_elem[irow * m_col + icol];
            mat[irow + m_row * (icol + 1)] = m_elem[irow * m_col + icol + 1];
            mat[irow + m_row * (icol + 2)] = m_elem[irow * m_col + icol + 2];
            mat[irow + m_row * (icol + 3)] = m_elem[irow * m_col + icol + 3];
        }

        for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
        {
            mat[irow + m_row * icol] = m_elem[irow * m_col + icol];
        }
    }

    return dtMatrix<m_col, m_row, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtMatrix<m_row, m_col, m_type>::Trace() const
{
    int cnt, i;
    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;
    m_type sum = 0;

    // for (uint16_t i = 0; i < num; i++)
    //     sum += m_elem[i * offset];
    for (cnt = num >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        sum += m_elem[i * offset];
        sum += m_elem[(i + 1) * offset];
        sum += m_elem[(i + 2) * offset];
        sum += m_elem[(i + 3) * offset];
    }
    for (cnt = num % 4; cnt > 0; cnt--, i++)
    {
        sum += m_elem[i * offset];
    }

    return sum;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtMatrix<m_row, m_col, m_type>::GetNorm() const
{
    uint16_t cnt, i = 0;
    m_type sqSum = 0;

    for (cnt = (m_row * m_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = (m_row * m_col) % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i];
    }

    return std::sqrt(sqSum);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtMatrix<m_row, m_col, m_type>::GetSqNorm() const
{
    uint16_t cnt, i = 0;
    m_type sqSum = 0;

    for (cnt = (m_row * m_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = (m_row * m_col) % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i];
    }

    return sqSum;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtMatrix<m_row, m_col, m_type>::Determinant() const
{
    return dtPartialPivLU<m_row, m_col, m_type>(*this).Determinant();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtNoPivLU<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::NoPivLU() const
{
    return dtNoPivLU<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtPartialPivLU<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::PartialPivLU() const
{
    return dtPartialPivLU<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLLT<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::LLT() const
{
    return dtLLT<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLDLT<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::LDLT() const
{
    return dtLDLT<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtQR<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::QR() const
{
    return dtQR<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtSVD<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::SVD() const
{
    return dtSVD<m_row, m_col, m_type>(*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::Inv(int8_t *isOk) const
{
    return dtPartialPivLU<m_row, m_col, m_type>(*this).Inverse(isOk);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_col, m_row, m_type> dtMatrix<m_row, m_col, m_type>::PInv(int8_t *isOk, m_type tolerance) const
{
    return dtSVD<m_row, m_col, m_type>(*this).Inverse(isOk, tolerance);
}

/* Assignment operators */
template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator=(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator+=(const dtMatrix<m_row, m_col, m_type> &m)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += m.m_elem[i];
        m_elem[i + 2] += m.m_elem[i + 2];
        m_elem[i + 4] += m.m_elem[i + 4];
        m_elem[i + 6] += m.m_elem[i + 6];
        m_elem[i + 1] += m.m_elem[i + 1];
        m_elem[i + 3] += m.m_elem[i + 3];
        m_elem[i + 5] += m.m_elem[i + 5];
        m_elem[i + 7] += m.m_elem[i + 7];
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += m.m_elem[i];
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator-=(const dtMatrix<m_row, m_col, m_type> &m)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= m.m_elem[i];
        m_elem[i + 2] -= m.m_elem[i + 2];
        m_elem[i + 4] -= m.m_elem[i + 4];
        m_elem[i + 6] -= m.m_elem[i + 6];
        m_elem[i + 1] -= m.m_elem[i + 1];
        m_elem[i + 3] -= m.m_elem[i + 3];
        m_elem[i + 5] -= m.m_elem[i + 5];
        m_elem[i + 7] -= m.m_elem[i + 7];
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= m.m_elem[i];
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator=(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator+=(const dtMatrix3<m_type, m_row, m_col> &m)
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator-=(const dtMatrix3<m_type, m_row, m_col> &m)
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator=(const dtRotation<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator+=(const dtRotation<m_type, m_row, m_col> &m)
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator-=(const dtRotation<m_type, m_row, m_col> &m)
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator=(const dtTransform<m_type, m_row, m_col> &m)
{
    m_elem[0] = m.m_R.m_elem[0];
    m_elem[1] = m.m_R.m_elem[1];
    m_elem[2] = m.m_R.m_elem[2];
    m_elem[3] = m.m_p.m_elem[0];
    m_elem[4] = m.m_R.m_elem[3];
    m_elem[5] = m.m_R.m_elem[4];
    m_elem[6] = m.m_R.m_elem[5];
    m_elem[7] = m.m_p.m_elem[1];
    m_elem[8] = m.m_R.m_elem[6];
    m_elem[9] = m.m_R.m_elem[7];
    m_elem[10] = m.m_R.m_elem[8];
    m_elem[11] = m.m_p.m_elem[2];
    m_elem[12] = 0;
    m_elem[13] = 0;
    m_elem[14] = 0;
    m_elem[15] = 1;

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator+=(const dtTransform<m_type, m_row, m_col> &m)
{
    m_elem[0] += m.m_R.m_elem[0];
    m_elem[1] += m.m_R.m_elem[1];
    m_elem[2] += m.m_R.m_elem[2];
    m_elem[3] += m.m_p.m_elem[0];
    m_elem[4] += m.m_R.m_elem[3];
    m_elem[5] += m.m_R.m_elem[4];
    m_elem[6] += m.m_R.m_elem[5];
    m_elem[7] += m.m_p.m_elem[1];
    m_elem[8] += m.m_R.m_elem[6];
    m_elem[9] += m.m_R.m_elem[7];
    m_elem[10] += m.m_R.m_elem[8];
    m_elem[11] += m.m_p.m_elem[2];
    m_elem[15] += 1;

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator-=(const dtTransform<m_type, m_row, m_col> &m)
{
    m_elem[0] -= m.m_R.m_elem[0];
    m_elem[1] -= m.m_R.m_elem[1];
    m_elem[2] -= m.m_R.m_elem[2];
    m_elem[3] -= m.m_p.m_elem[0];
    m_elem[4] -= m.m_R.m_elem[3];
    m_elem[5] -= m.m_R.m_elem[4];
    m_elem[6] -= m.m_R.m_elem[5];
    m_elem[7] -= m.m_p.m_elem[1];
    m_elem[8] -= m.m_R.m_elem[6];
    m_elem[9] -= m.m_R.m_elem[7];
    m_elem[10] -= m.m_R.m_elem[8];
    m_elem[11] -= m.m_p.m_elem[2];
    m_elem[15] -= 1;

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator=(const m_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] = s;
        m_elem[i + 2] = s;
        m_elem[i + 4] = s;
        m_elem[i + 6] = s;
        m_elem[i + 1] = s;
        m_elem[i + 3] = s;
        m_elem[i + 5] = s;
        m_elem[i + 7] = s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = s;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator+=(const m_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += s;
        m_elem[i + 2] += s;
        m_elem[i + 4] += s;
        m_elem[i + 6] += s;
        m_elem[i + 1] += s;
        m_elem[i + 3] += s;
        m_elem[i + 5] += s;
        m_elem[i + 7] += s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += s;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator-=(const m_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= s;
        m_elem[i + 2] -= s;
        m_elem[i + 4] -= s;
        m_elem[i + 6] -= s;
        m_elem[i + 1] -= s;
        m_elem[i + 3] -= s;
        m_elem[i + 5] -= s;
        m_elem[i + 7] -= s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= s;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator*=(const m_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] *= s;
        m_elem[i + 2] *= s;
        m_elem[i + 4] *= s;
        m_elem[i + 6] *= s;
        m_elem[i + 1] *= s;
        m_elem[i + 3] *= s;
        m_elem[i + 5] *= s;
        m_elem[i + 7] *= s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] *= s;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> &dtMatrix<m_row, m_col, m_type>::operator/=(const m_type s)
{
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] /= scalar;
        m_elem[i + 2] /= scalar;
        m_elem[i + 4] /= scalar;
        m_elem[i + 6] /= scalar;
        m_elem[i + 1] /= scalar;
        m_elem[i + 3] /= scalar;
        m_elem[i + 5] /= scalar;
        m_elem[i + 7] /= scalar;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] /= scalar;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCommaInit<m_row * m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row * m_col, m_type>(m_elem);
}

/* Arithmetic operators */
template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-() const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = -m_elem[i];
        mat[i + 2] = -m_elem[i + 2];
        mat[i + 4] = -m_elem[i + 4];
        mat[i + 6] = -m_elem[i + 6];
        mat[i + 1] = -m_elem[i + 1];
        mat[i + 3] = -m_elem[i + 3];
        mat[i + 5] = -m_elem[i + 5];
        mat[i + 7] = -m_elem[i + 7];
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = -m_elem[i];
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator+(const dtMatrix<m_row, m_col, m_type> &m) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] + m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] + m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] + m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] + m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] + m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] + m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] + m.m_elem[i + 7];
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-(const dtMatrix<m_row, m_col, m_type> &m) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] - m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] - m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] - m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] - m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] - m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] - m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] - m.m_elem[i + 7];
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator+(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] + m.m_elem[0];
    mat[1] = m_elem[1] + m.m_elem[1];
    mat[2] = m_elem[2] + m.m_elem[2];
    mat[3] = m_elem[3] + m.m_elem[3];
    mat[4] = m_elem[4] + m.m_elem[4];
    mat[5] = m_elem[5] + m.m_elem[5];
    mat[6] = m_elem[6] + m.m_elem[6];
    mat[7] = m_elem[7] + m.m_elem[7];
    mat[8] = m_elem[8] + m.m_elem[8];

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-(const dtMatrix3<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] - m.m_elem[0];
    mat[1] = m_elem[1] - m.m_elem[1];
    mat[2] = m_elem[2] - m.m_elem[2];
    mat[3] = m_elem[3] - m.m_elem[3];
    mat[4] = m_elem[4] - m.m_elem[4];
    mat[5] = m_elem[5] - m.m_elem[5];
    mat[6] = m_elem[6] - m.m_elem[6];
    mat[7] = m_elem[7] - m.m_elem[7];
    mat[8] = m_elem[8] - m.m_elem[8];

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator+(const dtRotation<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] + m.m_elem[0];
    mat[1] = m_elem[1] + m.m_elem[1];
    mat[2] = m_elem[2] + m.m_elem[2];
    mat[3] = m_elem[3] + m.m_elem[3];
    mat[4] = m_elem[4] + m.m_elem[4];
    mat[5] = m_elem[5] + m.m_elem[5];
    mat[6] = m_elem[6] + m.m_elem[6];
    mat[7] = m_elem[7] + m.m_elem[7];
    mat[8] = m_elem[8] + m.m_elem[8];

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-(const dtRotation<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] - m.m_elem[0];
    mat[1] = m_elem[1] - m.m_elem[1];
    mat[2] = m_elem[2] - m.m_elem[2];
    mat[3] = m_elem[3] - m.m_elem[3];
    mat[4] = m_elem[4] - m.m_elem[4];
    mat[5] = m_elem[5] - m.m_elem[5];
    mat[6] = m_elem[6] - m.m_elem[6];
    mat[7] = m_elem[7] - m.m_elem[7];
    mat[8] = m_elem[8] - m.m_elem[8];

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator+(const dtTransform<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] + m.m_R.m_elem[0];
    mat[1] = m_elem[1] + m.m_R.m_elem[1];
    mat[2] = m_elem[2] + m.m_R.m_elem[2];
    mat[3] = m_elem[3] + m.m_p.m_elem[0];
    mat[4] = m_elem[4] + m.m_R.m_elem[3];
    mat[5] = m_elem[5] + m.m_R.m_elem[4];
    mat[6] = m_elem[6] + m.m_R.m_elem[5];
    mat[7] = m_elem[7] + m.m_p.m_elem[1];
    mat[8] = m_elem[8] + m.m_R.m_elem[6];
    mat[9] = m_elem[9] + m.m_R.m_elem[7];
    mat[10] = m_elem[10] + m.m_R.m_elem[8];
    mat[11] = m_elem[11] + m.m_p.m_elem[2];
    mat[12] = m_elem[12];
    mat[13] = m_elem[13];
    mat[14] = m_elem[14];
    mat[15] = m_elem[15] + 1;

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-(const dtTransform<m_type, m_row, m_col> &m) const
{
    m_type mat[m_row * m_col];

    mat[0] = m_elem[0] - m.m_R.m_elem[0];
    mat[1] = m_elem[1] - m.m_R.m_elem[1];
    mat[2] = m_elem[2] - m.m_R.m_elem[2];
    mat[3] = m_elem[3] - m.m_p.m_elem[0];
    mat[4] = m_elem[4] - m.m_R.m_elem[3];
    mat[5] = m_elem[5] - m.m_R.m_elem[4];
    mat[6] = m_elem[6] - m.m_R.m_elem[5];
    mat[7] = m_elem[7] - m.m_p.m_elem[1];
    mat[8] = m_elem[8] - m.m_R.m_elem[6];
    mat[9] = m_elem[9] - m.m_R.m_elem[7];
    mat[10] = m_elem[10] - m.m_R.m_elem[8];
    mat[11] = m_elem[11] - m.m_p.m_elem[2];
    mat[12] = m_elem[12];
    mat[13] = m_elem[13];
    mat[14] = m_elem[14];
    mat[15] = m_elem[15] - 1;

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator+(const m_type s) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] + s;
        mat[i + 2] = m_elem[i + 2] + s;
        mat[i + 4] = m_elem[i + 4] + s;
        mat[i + 6] = m_elem[i + 6] + s;
        mat[i + 1] = m_elem[i + 1] + s;
        mat[i + 3] = m_elem[i + 3] + s;
        mat[i + 5] = m_elem[i + 5] + s;
        mat[i + 7] = m_elem[i + 7] + s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] + s;
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator-(const m_type s) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] - s;
        mat[i + 2] = m_elem[i + 2] - s;
        mat[i + 4] = m_elem[i + 4] - s;
        mat[i + 6] = m_elem[i + 6] - s;
        mat[i + 1] = m_elem[i + 1] - s;
        mat[i + 3] = m_elem[i + 3] - s;
        mat[i + 5] = m_elem[i + 5] - s;
        mat[i + 7] = m_elem[i + 7] - s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] - s;
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const m_type s) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] * s;
        mat[i + 2] = m_elem[i + 2] * s;
        mat[i + 4] = m_elem[i + 4] * s;
        mat[i + 6] = m_elem[i + 6] * s;
        mat[i + 1] = m_elem[i + 1] * s;
        mat[i + 3] = m_elem[i + 3] * s;
        mat[i + 5] = m_elem[i + 5] * s;
        mat[i + 7] = m_elem[i + 7] * s;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] * s;
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator/(const m_type s) const
{
    m_type mat[m_row * m_col];
    uint16_t cnt, i = 0;
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] / scalar;
        mat[i + 2] = m_elem[i + 2] / scalar;
        mat[i + 4] = m_elem[i + 4] / scalar;
        mat[i + 6] = m_elem[i + 6] / scalar;
        mat[i + 1] = m_elem[i + 1] / scalar;
        mat[i + 3] = m_elem[i + 3] / scalar;
        mat[i + 5] = m_elem[i + 5] / scalar;
        mat[i + 7] = m_elem[i + 7] / scalar;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] / scalar;
    }

    return dtMatrix(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtMatrix<m_col, col, m_type> &m) const
{
    m_type mat[m_row * col] = {
        0,
    };
    uint16_t cnt;
    uint16_t irow, icol, i;

    for (irow = 0; irow < m_row; ++irow)
    {
        for (icol = 0; icol < col; ++icol)
        {
            for (cnt = m_col >> 2u, i = 0; cnt > 0u; cnt--, i += 4u)
            {
                mat[irow * col + icol] += m_elem[irow * m_col + i] * m.m_elem[i * col + icol];
                mat[irow * col + icol] += m_elem[irow * m_col + i + 1] * m.m_elem[(i + 1) * col + icol];
                mat[irow * col + icol] += m_elem[irow * m_col + i + 2] * m.m_elem[(i + 2) * col + icol];
                mat[irow * col + icol] += m_elem[irow * m_col + i + 3] * m.m_elem[(i + 3) * col + icol];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, i++)
            {
                mat[irow * col + icol] += m_elem[irow * m_col + i] * m.m_elem[i * col + icol];
            }
        }
    }

    return dtMatrix<m_row, col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtMatrix3<m_type, m_col, m_col> &m) const
{
    m_type mat[m_row * m_col] = {
        0,
    };

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        for (uint16_t icol = 0; icol < m_col; ++icol)
        {
            mat[irow * m_col + icol] += m_elem[irow * m_col] * m.m_elem[icol];
            mat[irow * m_col + icol] += m_elem[irow * m_col + 1] * m.m_elem[3 + icol];
            mat[irow * m_col + icol] += m_elem[irow * m_col + 2] * m.m_elem[6 + icol];
        }
    }

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtRotation<m_type, m_col, m_col> &m) const
{
    m_type mat[m_row * m_col] = {
        0,
    };

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        for (uint16_t icol = 0; icol < m_col; ++icol)
        {
            mat[irow * m_col + icol] += m_elem[irow * m_col] * m.m_elem[icol];
            mat[irow * m_col + icol] += m_elem[irow * m_col + 1] * m.m_elem[3 + icol];
            mat[irow * m_col + icol] += m_elem[irow * m_col + 2] * m.m_elem[6 + icol];
        }
    }

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtTransform<m_type, m_col, m_col> &m) const
{
    m_type mat[m_row * m_col] = {
        0,
    };

    for (uint16_t irow = 0; irow < m_row * 4; irow += 4)
    {
        for (uint16_t icol = 0; icol < 3; ++icol)
        {
            mat[irow + icol] =
                m_elem[irow] * m.m_R.m_elem[icol] +
                m_elem[irow + 1] * m.m_R.m_elem[icol + 3] +
                m_elem[irow + 2] * m.m_R.m_elem[icol + 6];
        }
        mat[irow + 3] =
            m_elem[irow] * m.m_p.m_elem[0] +
            m_elem[irow + 1] * m.m_p.m_elem[1] +
            m_elem[irow + 2] * m.m_p.m_elem[2] +
            m_elem[irow + 3];
    }

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtVector<m_col, m_type> &v) const
{
    m_type vec[m_row] = {
        0,
    };
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < m_row; ++irow)
    {
        for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            vec[irow] += m_elem[irow * m_col + icol] * v.m_elem[icol];
            vec[irow] += m_elem[irow * m_col + icol + 1] * v.m_elem[icol + 1];
            vec[irow] += m_elem[irow * m_col + icol + 2] * v.m_elem[icol + 2];
            vec[irow] += m_elem[irow * m_col + icol + 3] * v.m_elem[icol + 3];
        }

        for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            vec[irow] += m_elem[irow * m_col + icol] * v.m_elem[icol];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtVector3<m_type, m_col> &v) const
{
    m_type vec[m_row];

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec[irow] = m_elem[irow * 3] * v.m_elem[0];
        vec[irow] += m_elem[irow * 3 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 3 + 2] * v.m_elem[2];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtVector4<m_type, m_col> &v) const
{
    m_type vec[m_row];

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec[irow] = m_elem[irow * 4] * v.m_elem[0];
        vec[irow] += m_elem[irow * 4 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 4 + 2] * v.m_elem[2];
        vec[irow] += m_elem[irow * 4 + 3] * v.m_elem[3];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtMatrix<m_row, m_col, m_type>::operator*(const dtVector6<m_type, m_col> &v) const
{
    m_type vec[m_row];

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec[irow] = m_elem[irow * 6] * v.m_elem[0];
        vec[irow] += m_elem[irow * 6 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 6 + 2] * v.m_elem[2];
        vec[irow] += m_elem[irow * 6 + 3] * v.m_elem[3];
        vec[irow] += m_elem[irow * 6 + 4] * v.m_elem[4];
        vec[irow] += m_elem[irow * 6 + 5] * v.m_elem[5];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator&(const dtVector<m_col, m_type> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    static_assert(m_row == 3, "This method is only for 3 x 3 matrix");
    static_assert(m_col == 3, "This method is only for 3 x 3 matrix");

    m_type mat[m_row * m_col];

    mat[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtMatrix<m_row, m_col, m_type>::operator&(const dtVector3<m_type, m_col> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    static_assert(m_row == 3, "This method is only for 3 x 3 matrix");
    static_assert(m_col == 3, "This method is only for 3 x 3 matrix");

    m_type mat[m_row * m_col];

    mat[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return dtMatrix<m_row, m_col, m_type>(mat);
}

/* Comparison operators */
template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator==(const dtMatrix<m_row, m_col, m_type> &m) const
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance)
            return false;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance)
            return false;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance)
            return false;
    }

    return true;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator!=(const dtMatrix<m_row, m_col, m_type> &m) const
{
    uint16_t cnt, i = 0;

    for (cnt = (m_row * m_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance)
            return true;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance)
            return true;
    }

    for (cnt = (m_row * m_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance)
            return true;
    }

    return false;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator==(const dtMatrix3<m_type, m_row, m_col> &m) const
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator!=(const dtMatrix3<m_type, m_row, m_col> &m) const
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator==(const dtRotation<m_type, m_row, m_col> &m) const
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator!=(const dtRotation<m_type, m_row, m_col> &m) const
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

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator==(const dtTransform<m_type, m_row, m_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance)
        return false;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance)
        return false;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance)
        return false;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance)
        return false;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance)
        return false;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance)
        return false;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance)
        return false;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance)
        return false;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance)
        return false;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance)
        return false;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance)
        return false;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance)
        return false;

    if (std::abs(m_elem[12]) > m_tolerance)
        return false;
    if (std::abs(m_elem[13]) > m_tolerance)
        return false;
    if (std::abs(m_elem[14]) > m_tolerance)
        return false;
    if (std::abs(1 - m_elem[15]) > m_tolerance)
        return false;

    return true;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline bool dtMatrix<m_row, m_col, m_type>::operator!=(const dtTransform<m_type, m_row, m_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance)
        return true;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance)
        return true;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance)
        return true;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance)
        return true;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance)
        return true;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance)
        return true;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance)
        return true;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance)
        return true;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance)
        return true;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance)
        return true;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance)
        return true;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance)
        return true;

    if (std::abs(m_elem[12]) > m_tolerance)
        return true;
    if (std::abs(m_elem[13]) > m_tolerance)
        return true;
    if (std::abs(m_elem[14]) > m_tolerance)
        return true;
    if (std::abs(1 - m_elem[15]) > m_tolerance)
        return true;

    return false;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtMatrix<m_row, m_col, m_type>::Print(const char endChar)
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
            printf("%7.6f ", (m_type)(m_elem[irow * m_col + icol]));
        }
        printf("\n");
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <uint16_t row, uint16_t col, typename type>
inline dtMatrix<row, col, type> operator*(const type s, const dtMatrix<row, col, type> &m)
{
    type mat[row * col];
    uint16_t cnt, i = 0;

    for (cnt = (row * col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m.m_elem[i] * s;
        mat[i + 2] = m.m_elem[i + 2] * s;
        mat[i + 4] = m.m_elem[i + 4] * s;
        mat[i + 6] = m.m_elem[i + 6] * s;
        mat[i + 1] = m.m_elem[i + 1] * s;
        mat[i + 3] = m.m_elem[i + 3] * s;
        mat[i + 5] = m.m_elem[i + 5] * s;
        mat[i + 7] = m.m_elem[i + 7] * s;
    }

    for (cnt = (row * col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m.m_elem[i] * s;
    }

    return dtMatrix<row, col, type>(mat);
}

} // namespace dtMath

#endif // DTMATH_DTMATRIX_TPP_