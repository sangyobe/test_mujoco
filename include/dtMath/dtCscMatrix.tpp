/*!
\file       dtCscMatrix.h
\brief      dtMath, Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX_TPP_
#define DTMATH_DTCSC_MATRIX_TPP_

#include "dtCscMatrix.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type>::dtCscMatrix(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr)
{
    m_elemNum = elemNum;
    memcpy(m_elem, element, sizeof(m_type) * m_elemNum);
    memcpy(m_rowIdx, rowIdx, sizeof(m_rowIdx));
    memcpy(m_colPtr, colPtr, sizeof(m_colPtr));
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type>::dtCscMatrix()
{
    m_elemNum = 0;
    memset(m_elem, 0, sizeof(m_elem));
    memset(m_rowIdx, 0, sizeof(m_rowIdx));
    memset(m_colPtr, 0, sizeof(m_colPtr));
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type>::dtCscMatrix(const dtMatrix<m_row, m_col, m_type> &m)
{
    uint16_t icol;
    uint16_t i = 0;
    m_type tolerance = std::numeric_limits<m_type>::epsilon();

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(m.m_elem[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type>::dtCscMatrix(const dtCscMatrix &m)
{
    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(m_elem));
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(m_rowIdx));
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline const m_type *const dtCscMatrix<m_row, m_col, m_type>::GetDataAddr() const
{
    return m_elem;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline const int *const dtCscMatrix<m_row, m_col, m_type>::GetRowIdx() const
{
    return m_rowIdx;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline const int *const dtCscMatrix<m_row, m_col, m_type>::GetColPtr() const
{
    return m_colPtr;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtCscMatrix<m_row, m_col, m_type>::GetRowVec(const uint16_t idxRow) const
{
    m_type vec[m_col] = {
        0,
    };

    if (m_elemNum == 0)
        return dtVector<m_col, m_type>();

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
                vec[j] = m_elem[i];
        }
    }

    return dtVector<m_col, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtCscMatrix<m_row, m_col, m_type>::GetColVec(const uint16_t idxCol) const
{
    m_type vec[m_row] = {
        0,
    };

    if (m_elemNum == 0)
        return dtVector<m_row, m_type>();

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        vec[m_rowIdx[i]] = m_elem[i];
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtCscMatrix<m_row, m_col, m_type>::GetRowVec(const uint16_t idxRow, dtVector<m_col, m_type> &v) const
{
    int8_t rtn = -1;

    memset(v.m_elem, 0, sizeof(m_type) * m_col);

    if (m_elemNum == 0)
        return -1;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
            {
                v.m_elem[j] = m_elem[i];
                rtn = 0;
                break;
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtCscMatrix<m_row, m_col, m_type>::GetColVec(const uint16_t idxCol, dtVector<m_row, m_type> &v) const
{
    uint16_t i;

    memset(v.m_elem, 0, sizeof(m_type) * m_row);

    if (m_elemNum == 0)
        return -1;

    for (i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        v.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    if (m_colPtr[idxCol] == i)
        return -1;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtCscMatrix<m_row, m_col, m_type>::GetDenseMat() const
{
    m_type mat[m_row * m_col];

    memset(mat, 0, sizeof(mat));

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            mat[m_rowIdx[i] * m_col + j] = m_elem[i];
        }
    }

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_col, m_row, m_type> dtCscMatrix<m_row, m_col, m_type>::Transpose() const
{
    // Transpose ==  csc -> csr

    m_type elem[m_row * m_col];
    int colIdx[m_row * m_col] = {
        0,
    };
    int rowPtr[m_row + 1] = {
        0,
    }; // m_row + 1
    int tmpPtr[m_row + 1] = {
        0,
    };

    if (!m_elemNum)
        return dtMatrix<m_col, m_row, m_type>();

    // compute number of non-zero entries per row
    for (uint16_t n = 0; n < m_elemNum; n++)
    {
        rowPtr[m_rowIdx[n]]++;
    }

    // cumsum the elemNum per row to get rowPtr[]
    for (uint16_t i = 0, cumsum = 0, temp; i < m_row; i++)
    {
        temp = rowPtr[i]; // number of non-zero entries per row
        rowPtr[i] = cumsum;
        cumsum += temp;
    }
    rowPtr[m_row] = m_elemNum;
    memcpy(tmpPtr, rowPtr, sizeof(rowPtr));

    // compute column index and data element
    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j], ii, jj; i < m_colPtr[j + 1]; i++)
        {
            ii = m_rowIdx[i];
            jj = tmpPtr[ii];

            colIdx[jj] = j;
            elem[jj] = m_elem[i];

            tmpPtr[ii]++;
        }
    }

    return dtCscMatrix<m_col, m_row, m_type>(elem, m_elemNum, colIdx, rowPtr);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtCscMatrix<m_row, m_col, m_type>::GetNorm() const
{
    m_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i];

    return std::sqrt(sqSum);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtCscMatrix<m_row, m_col, m_type>::GetSqNorm() const
{
    m_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i];

    return sqSum;
}

/* Assignment operators */
template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> &dtCscMatrix<m_row, m_col, m_type>::operator=(const dtCscMatrix &m)
{
    m_elemNum = m.m_elem;
    memcpy(m_elem, m.m_elem, sizeof(m_elem));
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(m_rowIdx));
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> &dtCscMatrix<m_row, m_col, m_type>::operator*=(const m_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] *= s;
        m_elem[i + 2] *= s;
        m_elem[i + 1] *= s;
        m_elem[i + 3] *= s;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] *= s;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> &dtCscMatrix<m_row, m_col, m_type>::operator/=(const m_type s)
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

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] /= scalar;
        m_elem[i + 2] /= scalar;
        m_elem[i + 1] /= scalar;
        m_elem[i + 3] /= scalar;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] /= scalar;
    }

    return (*this);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> dtCscMatrix<m_row, m_col, m_type>::operator*(const m_type s) const
{
    m_type elem[m_elemNum];
    uint16_t cnt, i = 0;

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        elem[i] = m_elem[i] * s;
        elem[i + 2] = m_elem[i + 2] * s;
        elem[i + 1] = m_elem[i + 1] * s;
        elem[i + 3] = m_elem[i + 3] * s;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        elem[i] = m_elem[i] * s;
    }

    return dtCscMatrix(elem, m_elemNum, m_rowIdx, m_colPtr);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtCscMatrix<m_row, m_col, m_type> dtCscMatrix<m_row, m_col, m_type>::operator/(const m_type s) const
{
    m_type elem[m_elemNum];
    uint16_t cnt, i = 0;
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        elem[i] = m_elem[i] / scalar;
        elem[i + 2] = m_elem[i + 2] / scalar;
        elem[i + 1] = m_elem[i + 1] / scalar;
        elem[i + 3] = m_elem[i + 3] / scalar;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        elem[i] = m_elem[i] / scalar;
    }

    return dtCscMatrix(elem, m_elemNum, m_rowIdx, m_colPtr);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_row, m_type> dtCscMatrix<m_row, m_col, m_type>::operator*(const dtVector<m_col, m_type> &v) const
{
    m_type vec[m_row] = {
        0,
    };

    if (!m_elemNum)
        return dtVector<m_row, m_type>();

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return dtVector<m_row, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtCscMatrix<m_row, m_col, m_type>::TposeVec(const dtVector<m_row, m_type> &v) const
{
    m_type vec[m_col] = {
        0,
    };

    if (!m_elemNum)
        return dtVector<m_col, m_type>();

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return dtVector<m_col, m_type>(vec);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtCscMatrix<m_row, m_col, m_type>::Print(const char endChar)
{
    printf("dat vec =");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%7.2f ", (m_type)(m_elem[i]));
    }
    printf("\nvec num = %d\n", m_elemNum);

    printf("row idx = ");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%d ", m_rowIdx[i]);
    }
    printf("\n");

    printf("col ptr = ");
    for (uint16_t i = 0; i < m_col + 1; i++)
    {
        printf("%d ", m_colPtr[i]);
    }
    printf("\n%c", endChar);
}

} // namespace dtMath

#endif // DTMATH_DTCSC_MATRIX_TPP_
