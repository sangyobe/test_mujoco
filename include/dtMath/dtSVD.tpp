/*!
\file       dtSVD.h
\brief      dtMath, Singular Value Decomposition solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD_TPP_
#define DTMATH_DTSVD_TPP_

#include "dtSVD.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtSVD<m_row, m_col, m_type>::dtSVD()
{
    memset(m_A, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtSVD<m_row, m_col, m_type>::dtSVD(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
    {
        m_isOk = 0;
        return;
    }

    if (m_row >= m_col)
    {
        memcpy(m_A, element, n_byte);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn();
        m_isOk = 1;
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < m_row; ++irow)
        //     for (uint16_t icol = 0; icol < m_col; ++icol)
        //         m_A[icol * m_row + irow] = element[irow * m_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = element[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = element[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = element[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN();
        m_isOk = 1;
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtSVD<m_row, m_col, m_type>::dtSVD(const dtMatrix<m_row, m_col, m_type> &m)
{
    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(m_type) * m_row * m_col);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn();
        m_isOk = 1;
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < m_row; ++irow)
        //     for (uint16_t icol = 0; icol < m_col; ++icol)
        //         m_A[icol * m_row + irow] = m.m_elem[irow * m_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN();
        m_isOk = 1;
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtSVD<m_row, m_col, m_type>::dtSVD(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_A, m.m_elem, sizeof(m_type) * m_row * m_col);

    // Compute SVD
    HouseholdersReductionToBidiagonalForm_Mxn();
    if (GivensReductionToDiagonalForm_Mxn())
    {
        m_isOk = 0;
        return;
    }
    SortByDecreasingSingularValues_Mxn();
    m_isOk = 1;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    if (m_row >= m_col)
    {
        memcpy(m_A, element, n_byte);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn();
        m_isOk = 1;
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < m_row; ++irow)
        //     for (uint16_t icol = 0; icol < m_col; ++icol)
        //         m_A[icol * m_row + irow] = element[irow * m_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = element[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = element[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = element[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN();
        m_isOk = 1;
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(m_type) * m_row * m_col);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn();
        m_isOk = 1;
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < m_row; ++irow)
        //     for (uint16_t icol = 0; icol < m_col; ++icol)
        //         m_A[icol * m_row + irow] = m.m_elem[irow * m_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN();
        m_isOk = 1;
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_A, m.m_elem, sizeof(m_type) * m_row * m_col);

    // Compute SVD
    HouseholdersReductionToBidiagonalForm_Mxn();
    if (GivensReductionToDiagonalForm_Mxn())
    {
        m_isOk = 0;
        return -1;
    }
    SortByDecreasingSingularValues_Mxn();
    m_isOk = 1;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_row, m_type> dtSVD<m_row, m_col, m_type>::GetMatrixU() const
{
    if (m_row > m_col)
    {
        m_type U[m_row * m_row] = {
            0,
        };

        for (int i = 0; i < m_row; i++)
            memcpy(&U[i * m_row], &m_U[i * m_col], sizeof(m_type) * m_col);

        return dtMatrix<m_row, m_row, m_type>(U);
    }
    else
    {
        return dtMatrix<m_row, m_row, m_type>(m_U);
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtSVD<m_row, m_col, m_type>::GetMatrixS() const
{
    return dtMatrix<m_row, m_col, m_type>(m_S, 'd');
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_col, m_col, m_type> dtSVD<m_row, m_col, m_type>::GetMatrixV() const
{
    if (m_row < m_col)
    {
        m_type V[m_col * m_col] = {
            0,
        };

        for (int i = 0; i < m_col; i++)
            memcpy(&V[i * m_col], &m_V[i * m_row], sizeof(m_type) * m_row);

        return dtMatrix<m_col, m_col, m_type>(V);
    }
    else
    {
        return dtMatrix<m_col, m_col, m_type>(m_V);
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline int8_t dtSVD<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, dtMatrix<m_col, col, m_type> &x, m_type tolerance)
{
    int i, j, k, c;
    m_type *pU, *pV;
    m_type sum;

    if (!m_isOk)
        return -1;

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV;
    m_type sum;

    if (!m_isOk)
        return -1;

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline dtMatrix<m_col, col, m_type> dtSVD<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, int8_t *isOk, m_type tolerance)
{
    int i, j, k, c;
    m_type *pU, *pV;
    m_type sum;
    m_type x[m_col * col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtMatrix<m_col, col, m_type>();
    }

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
            {
                x[i * col + c] = 0;
                for (j = 0; j < m_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
            {
                x[i * col + c] = 0;
                for (j = 0; j < m_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    return dtMatrix<m_col, col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtSVD<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, int8_t *isOk, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV;
    m_type sum;
    m_type x[m_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
        {
            x[i] = 0;
            for (j = 0; j < m_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
        {
            x[i] = 0;
            for (j = 0; j < m_row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                        sum += *(pU + j) * b.m_elem[k];

                    x[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Inverse(dtMatrix<m_col, m_row, m_type> &inv, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV, *pInvA;
    m_type tol;

    if (!m_isOk)
        return -1;

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::Inverse(dtMatrix3<m_type, m_col, m_row> &inv, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV, *pInvA;
    m_type tol;

    if (!m_isOk)
        return -1;

    tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
    if (tolerance < tol)
        tolerance = tol;

    for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_col)
    {
        for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
        {
            for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                if (m_S[k] > tolerance)
                    *pInvA += *(pV + k) * *pU / m_S[k];
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_col, m_row, m_type> dtSVD<m_row, m_col, m_type>::Inverse(int8_t *isOk, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV, *pInvA;
    m_type tol;

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtMatrix<m_col, m_row, m_type>();
    }

    memset(m_inv, 0, sizeof(m_type) * m_row * m_col);

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return dtMatrix<m_col, m_row, m_type>(m_inv);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::InverseArray(m_type *inv, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV, *pInvA;
    m_type tol;

    if (!m_isOk)
        return -1;

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type *dtSVD<m_row, m_col, m_type>::InverseArray(int8_t *isOk, m_type tolerance)
{
    int i, j, k;
    m_type *pU, *pV, *pInvA;
    m_type tol;
    memset(m_inv, 0, sizeof(m_type) * m_row * m_col);

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return m_inv;
    }

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<m_type>::epsilon() * m_S[0] * (m_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return m_inv;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtSVD<m_row, m_col, m_type>::HouseholdersReductionToBidiagonalForm_Mxn()
{
    int i, j, k, ip1;
    m_type s, s2, si, scale;
    m_type *pU, *pUi, *pV, *pVi;
    m_type halfNormSquared;

    // Copy A to U
    memcpy(m_U, m_A, sizeof(m_type) * m_row * m_col);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pUi = m_U, ip1 = 1; i < m_col; pUi += m_col, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pU = pUi, scale = 0; j < m_row; j++, pU += m_col)
            scale += std::abs(*(pU + i));

        if (scale > std::numeric_limits<m_type>::epsilon())
        {
            for (j = i, pU = pUi, s2 = 0; j < m_row; j++, pU += m_col)
            {
                *(pU + i) /= scale;
                s2 += *(pU + i) * *(pU + i);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pUi + i) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/u'u
            halfNormSquared = *(pUi + i) * s - s2;

            // Transform remaining columns by the Householder transform.
            *(pUi + i) -= s;

            for (j = ip1; j < m_col; j++)
            {
                for (k = i, si = 0, pU = pUi; k < m_row; k++, pU += m_col)
                    si += *(pU + i) * *(pU + j);

                si /= halfNormSquared;

                for (k = i, pU = pUi; k < m_row; k++, pU += m_col)
                    *(pU + j) += si * *(pU + i);
            }
        }

        for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
            *(pU + i) *= scale;

        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= m_row || i == (m_col - 1))
            continue;

        for (j = ip1; j < m_col; j++)
            scale += std::abs(*(pUi + j));

        if (scale > std::numeric_limits<m_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < m_col; j++)
            {
                *(pUi + j) /= scale;
                s2 += *(pUi + j) * *(pUi + j);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pUi + ip1) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/u'u
            halfNormSquared = *(pUi + ip1) * s - s2;

            // Transform the rows by the Householder transform.
            *(pUi + ip1) -= s;

            for (k = ip1; k < m_col; k++)
                m_superDiagonal[k] = *(pUi + k) / halfNormSquared;

            if (i < (m_row - 1))
            {
                for (j = ip1, pU = pUi + m_col; j < m_row; j++, pU += m_col)
                {
                    for (k = ip1, si = 0; k < m_col; k++)
                        si += *(pUi + k) * *(pU + k);

                    for (k = ip1; k < m_col; k++)
                        *(pU + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < m_col; k++)
                *(pUi + k) *= scale;
        }
    }

    /* Update V */
    pUi = m_U + m_col * (m_col - 2);
    pVi = m_V + m_col * (m_col - 1);
    *(pVi + m_col - 1) = 1;
    s = m_superDiagonal[m_col - 1];
    pVi -= m_col;

    for (i = m_col - 2, ip1 = m_col - 1; i >= 0; i--, pUi -= m_col, pVi -= m_col, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<m_type>::epsilon())
        {
            pV = pVi + m_col;

            for (j = ip1; j < m_col; j++, pV += m_col)
                *(pV + i) = (*(pUi + j) / *(pUi + ip1)) / s;

            for (j = ip1; j < m_col; j++)
            {
                si = 0;

                for (k = ip1, pV = pVi + m_col; k < m_col; k++, pV += m_col)
                    si += *(pUi + k) * *(pV + j);

                for (k = ip1, pV = pVi + m_col; k < m_col; k++, pV += m_col)
                    *(pV + j) += si * *(pV + i);
            }
        }

        pV = pVi + m_col;
        for (j = ip1; j < m_col; j++, pV += m_col)
        {
            *(pVi + j) = 0;
            *(pV + i) = 0;
        }

        *(pVi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update U */
    pUi = m_U + m_col * (m_col - 1);
    for (i = m_col - 1, ip1 = m_col; i >= 0; ip1 = i, i--, pUi -= m_col)
    {
        s = m_S[i];

        for (j = ip1; j < m_col; j++)
            *(pUi + j) = 0;

        if (std::abs(s) > std::numeric_limits<m_type>::epsilon())
        {
            for (j = ip1; j < m_col; j++)
            {
                si = 0;
                pU = pUi + m_col;

                for (k = ip1; k < m_row; k++, pU += m_col)
                    si += *(pU + i) * *(pU + j);

                si = (si / *(pUi + i)) / s;
                for (k = i, pU = pUi; k < m_row; k++, pU += m_col)
                    *(pU + j) += si * *(pU + i);
            }

            for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
                *(pU + i) /= s;
        }
        else
        {
            for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
                *(pU + i) = 0;
        }

        *(pUi + i) += 1;
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::GivensReductionToDiagonalForm_Mxn()
{
    m_type epsilon;
    m_type c, s;
    m_type f, g, h;
    m_type x, y, z;
    m_type *pU, *pV;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < m_col; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<m_type>::epsilon();

    for (k = m_col - 1; k >= 0; k--)
    {
        iteration_count = 0;
        while (1)
        {
            rotation_test = 1;
            for (m = k; m >= 0; m--)
            {
                if (std::abs(m_superDiagonal[m]) <= epsilon)
                {
                    rotation_test = 0;
                    break;
                }

                if (std::abs(m_S[m - 1]) <= epsilon)
                    break;
            }

            if (rotation_test)
            {
                c = 0;
                s = 1;
                for (i = m; i <= k; i++)
                {
                    f = s * m_superDiagonal[i];
                    m_superDiagonal[i] *= c;

                    if (std::abs(f) <= epsilon)
                        break;

                    g = m_S[i];
                    h = std::sqrt(f * f + g * g);
                    m_S[i] = h;
                    c = g / h;
                    s = -f / h;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_col)
                    {
                        y = *(pU + m - 1);
                        z = *(pU + i);
                        *(pU + m - 1) = y * c + z * s;
                        *(pU + i) = -y * s + z * c;
                    }
                }
            }

            z = m_S[k];

            if (m == k)
            {
                if (z < 0)
                {
                    m_S[k] = -z;
                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_col)
                        *(pV + k) = -*(pV + k);
                }
                break;
            }
            else
            {
                if (iteration_count >= MAX_ITERATION_COUNT)
                    return -1;

                iteration_count++;
                x = m_S[m];
                y = m_S[k - 1];
                g = m_superDiagonal[k - 1];
                h = m_superDiagonal[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                g = std::sqrt(f * f + 1);

                if (f < 0)
                    g = -g;

                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                // Next QR Transformtion
                c = 1;
                s = 1;

                for (i = m + 1; i <= k; i++)
                {
                    g = m_superDiagonal[i];
                    y = m_S[i];
                    h = s * g;
                    g *= c;
                    z = std::sqrt(f * f + h * h);
                    m_superDiagonal[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y *= c;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_col)
                    {
                        x = *(pV + i - 1);
                        z = *(pV + i);
                        *(pV + i - 1) = x * c + z * s;
                        *(pV + i) = -x * s + z * c;
                    }

                    z = std::sqrt(f * f + h * h);
                    m_S[i - 1] = z;

                    if (std::abs(z) > std::numeric_limits<m_type>::epsilon())
                    {
                        c = f / z;
                        s = h / z;
                    }

                    f = c * g + s * y;
                    x = -s * g + c * y;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_col)
                    {
                        y = *(pU + i - 1);
                        z = *(pU + i);
                        *(pU + i - 1) = c * y + s * z;
                        *(pU + i) = -s * y + c * z;
                    }
                }
                m_superDiagonal[m] = 0;
                m_superDiagonal[k] = f;
                m_S[k] = x;
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtSVD<m_row, m_col, m_type>::SortByDecreasingSingularValues_Mxn()
{
    int i, j, maxIdx;
    m_type temp;
    m_type *pM1, *pM2;

    for (i = 0; i < m_col - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < m_col; j++)
        {
            if (m_S[j] > m_S[maxIdx])
                maxIdx = j;
        }

        if (maxIdx == i)
            continue;

        temp = m_S[i];
        m_S[i] = m_S[maxIdx];
        m_S[maxIdx] = temp;
        pM1 = m_U + maxIdx;
        pM2 = m_U + i;

        for (j = 0; j < m_row; j++, pM1 += m_col, pM2 += m_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_V + maxIdx;
        pM2 = m_V + i;

        for (j = 0; j < m_col; j++, pM1 += m_col, pM2 += m_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtSVD<m_row, m_col, m_type>::HouseholdersReductionToBidiagonalForm_mxN()
{
    int i, j, k, ip1;
    m_type s, s2, si, scale;
    m_type *pV, *pVi, *pU, *pUi;
    m_type halfNormSquared;

    // Copy A to V
    memcpy(m_V, m_A, sizeof(m_type) * m_col * m_row);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pVi = m_V, ip1 = 1; i < m_row; pVi += m_row, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pV = pVi, scale = 0; j < m_col; j++, pV += m_row)
            scale += std::abs(*(pV + i));

        if (scale > std::numeric_limits<m_type>::epsilon())
        {
            for (j = i, pV = pVi, s2 = 0; j < m_col; j++, pV += m_row)
            {
                *(pV + i) /= scale;
                s2 += *(pV + i) * *(pV + i);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pVi + i) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/v'v
            halfNormSquared = *(pVi + i) * s - s2;

            // Transform remaining columns by the Householder transform.
            *(pVi + i) -= s;

            for (j = ip1; j < m_row; j++)
            {
                for (k = i, si = 0, pV = pVi; k < m_col; k++, pV += m_row)
                    si += *(pV + i) * *(pV + j);

                si /= halfNormSquared;

                for (k = i, pV = pVi; k < m_col; k++, pV += m_row)
                    *(pV + j) += si * *(pV + i);
            }
        }
        for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
            *(pV + i) *= scale;
        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= m_col || i == (m_row - 1))
            continue;

        for (j = ip1; j < m_row; j++)
            scale += std::abs(*(pVi + j));

        if (scale > std::numeric_limits<m_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < m_row; j++)
            {
                *(pVi + j) /= scale;
                s2 += *(pVi + j) * *(pVi + j);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pVi + ip1) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/v'v
            halfNormSquared = *(pVi + ip1) * s - s2;

            // Transform the rows by the Householder transform.
            *(pVi + ip1) -= s;

            for (k = ip1; k < m_row; k++)
                m_superDiagonal[k] = *(pVi + k) / halfNormSquared;

            if (i < (m_col - 1))
            {
                for (j = ip1, pV = pVi + m_row; j < m_col; j++, pV += m_row)
                {
                    for (k = ip1, si = 0; k < m_row; k++)
                        si += *(pVi + k) * *(pV + k);

                    for (k = ip1; k < m_row; k++)
                        *(pV + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < m_row; k++)
                *(pVi + k) *= scale;
        }
    }

    /* Update U */
    pVi = m_V + m_row * (m_row - 2);
    pUi = m_U + m_row * (m_row - 1);
    *(pUi + m_row - 1) = 1;
    s = m_superDiagonal[m_row - 1];
    pUi -= m_row;

    for (i = m_row - 2, ip1 = m_row - 1; i >= 0; i--, pVi -= m_row, pUi -= m_row, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<m_type>::epsilon())
        {
            pU = pUi + m_row;

            for (j = ip1; j < m_row; j++, pU += m_row)
                *(pU + i) = (*(pVi + j) / *(pVi + ip1)) / s;

            for (j = ip1; j < m_row; j++)
            {
                si = 0;

                for (k = ip1, pU = pUi + m_row; k < m_row; k++, pU += m_row)
                    si += *(pVi + k) * *(pU + j);

                for (k = ip1, pU = pUi + m_row; k < m_row; k++, pU += m_row)
                    *(pU + j) += si * *(pU + i);
            }
        }

        pU = pUi + m_row;
        for (j = ip1; j < m_row; j++, pU += m_row)
        {
            *(pUi + j) = 0;
            *(pU + i) = 0;
        }

        *(pUi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update V */
    pVi = m_V + m_row * (m_row - 1);
    for (i = m_row - 1, ip1 = m_row; i >= 0; ip1 = i, i--, pVi -= m_row)
    {
        s = m_S[i];

        for (j = ip1; j < m_row; j++)
            *(pVi + j) = 0;

        if (std::abs(s) > std::numeric_limits<m_type>::epsilon())
        {
            for (j = ip1; j < m_row; j++)
            {
                si = 0;
                pV = pVi + m_row;

                for (k = ip1; k < m_col; k++, pV += m_row)
                    si += *(pV + i) * *(pV + j);

                si = (si / *(pVi + i)) / s;
                for (k = i, pV = pVi; k < m_col; k++, pV += m_row)
                    *(pV + j) += si * *(pV + i);
            }

            for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
                *(pV + i) /= s;
        }
        else
        {
            for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
                *(pV + i) = 0;
        }

        *(pVi + i) += 1;
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtSVD<m_row, m_col, m_type>::GivensReductionToDiagonalForm_mxN()
{
    m_type epsilon;
    m_type c, s;
    m_type f, g, h;
    m_type x, y, z;
    m_type *pV, *pU;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < m_row; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<m_type>::epsilon();

    for (k = m_row - 1; k >= 0; k--)
    {
        iteration_count = 0;
        while (1)
        {
            rotation_test = 1;
            for (m = k; m >= 0; m--)
            {
                if (std::abs(m_superDiagonal[m]) <= epsilon)
                {
                    rotation_test = 0;
                    break;
                }

                if (std::abs(m_S[m - 1]) <= epsilon)
                    break;
            }

            if (rotation_test)
            {
                c = 0;
                s = 1;
                for (i = m; i <= k; i++)
                {
                    f = s * m_superDiagonal[i];
                    m_superDiagonal[i] *= c;

                    if (std::abs(f) <= epsilon)
                        break;

                    g = m_S[i];
                    h = std::sqrt(f * f + g * g);
                    m_S[i] = h;
                    c = g / h;
                    s = -f / h;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_row)
                    {
                        y = *(pV + m - 1);
                        z = *(pV + i);
                        *(pV + m - 1) = y * c + z * s;
                        *(pV + i) = -y * s + z * c;
                    }
                }
            }

            z = m_S[k];

            if (m == k)
            {
                if (z < 0)
                {
                    m_S[k] = -z;
                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_row)
                        *(pU + k) = -*(pU + k);
                }
                break;
            }
            else
            {
                if (iteration_count >= MAX_ITERATION_COUNT)
                    return -1;

                iteration_count++;
                x = m_S[m];
                y = m_S[k - 1];
                g = m_superDiagonal[k - 1];
                h = m_superDiagonal[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                g = std::sqrt(f * f + 1);

                if (f < 0)
                    g = -g;

                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                // Next QR Transformtion
                c = 1;
                s = 1;

                for (i = m + 1; i <= k; i++)
                {
                    g = m_superDiagonal[i];
                    y = m_S[i];
                    h = s * g;
                    g *= c;
                    z = std::sqrt(f * f + h * h);
                    m_superDiagonal[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y *= c;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_row)
                    {
                        x = *(pU + i - 1);
                        z = *(pU + i);
                        *(pU + i - 1) = x * c + z * s;
                        *(pU + i) = -x * s + z * c;
                    }

                    z = std::sqrt(f * f + h * h);
                    m_S[i - 1] = z;

                    if (std::abs(z) > std::numeric_limits<m_type>::epsilon())
                    {
                        c = f / z;
                        s = h / z;
                    }

                    f = c * g + s * y;
                    x = -s * g + c * y;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_row)
                    {
                        y = *(pV + i - 1);
                        z = *(pV + i);
                        *(pV + i - 1) = c * y + s * z;
                        *(pV + i) = -s * y + c * z;
                    }
                }
                m_superDiagonal[m] = 0;
                m_superDiagonal[k] = f;
                m_S[k] = x;
            }
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void dtSVD<m_row, m_col, m_type>::SortByDecreasingSingularValues_mxN()
{
    int i, j, maxIdx;
    m_type temp;
    m_type *pM1, *pM2;

    for (i = 0; i < m_row - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < m_row; j++)
        {
            if (m_S[j] > m_S[maxIdx])
                maxIdx = j;
        }

        if (maxIdx == i)
            continue;

        temp = m_S[i];
        m_S[i] = m_S[maxIdx];
        m_S[maxIdx] = temp;
        pM1 = m_V + maxIdx;
        pM2 = m_V + i;

        for (j = 0; j < m_col; j++, pM1 += m_row, pM2 += m_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_U + maxIdx;
        pM2 = m_U + i;

        for (j = 0; j < m_row; j++, pM1 += m_row, pM2 += m_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

} // namespace dtMath

#endif // DTMATH_DTSVD_TPP_
