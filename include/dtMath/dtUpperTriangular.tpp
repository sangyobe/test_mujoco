/*!
\file       dtUpperTriangular.h
\brief      dtMath, Upper triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTUPPER_TRIANGULAR_TPP_
#define DTMATH_DTUPPER_TRIANGULAR_TPP_

#include "dtUpperTriangular.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtUpperTriangular<m_row, m_col, m_type>::Solve(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    if (m_row != m_col)
        return -1;

    int k, i;
    m_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = m_col - 1, pU += m_row * (m_col - 1); k >= 0; pU -= m_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<m_type>::epsilon())
            return -1; // The matrix U is singular

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < m_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtUpperTriangular<m_row, m_col, m_type>::Solve(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, int8_t *isOk)
{
    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    int k, i;
    m_type x[m_col] = {
        0,
    };
    m_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = m_col - 1, pU += m_row * (m_col - 1); k >= 0; pU -= m_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<m_type>::epsilon())
        {
            // The matrix U is singular
            if (isOk)
                *isOk = 0;
            return dtVector<m_col, m_type>();
        }

        x[k] = b.m_elem[k];

        for (i = k + 1; i < m_col; i++)
            x[k] -= x[i] * *(pU + i);

        x[k] /= *(pU + k);
    }

    if (isOk)
        *isOk = 1;

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtUpperTriangular<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> U, dtMatrix<m_row, m_col, m_type> &invU)
{
    int i, j, k;
    m_type *pMi, *pMk;
    m_type sum;

    if (m_row != m_col)
        return -1;

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<m_type>::epsilon())
            return -1;
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix U, for row i */
    for (i = m_row - 2, pMi = U.m_elem + m_col * (m_row - 2); i >= 0; pMi -= m_col, i--)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + m_col; k <= j; pMk += m_col, k++)
            {
                sum += *(pMi + k) * *(pMk + j);
            }
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invU = U;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtUpperTriangular<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> U, int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMk;
    m_type sum;

    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtMatrix<m_row, m_col, m_type>();
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<m_type>::epsilon())
        {
            if (isOk)
                *isOk = 0;
            return dtMatrix<m_row, m_col, m_type>();
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = m_row - 2, pMi = U.m_elem + m_col * (m_row - 2); i >= 0; pMi -= m_col, i--)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + m_col; k <= j; pMk += m_col, k++)
            {
                sum += *(pMi + k) * *(pMk + j);
            }
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    if (isOk)
        *isOk = 1;

    return U;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtUpperTriangular<m_row, m_col, m_type>::SolveUnit(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    if (m_row != m_col)
        return -1;

    int k, i;
    m_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[m_col - 1] = b.m_elem[m_row - 1];
    for (k = m_col - 2, pU += m_row * (m_col - 2); k >= 0; pU -= m_col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < m_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtUpperTriangular<m_row, m_col, m_type>::SolveUnit(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, int8_t *isOk)
{
    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    int k, i;
    m_type x[m_col] = {
        0,
    };
    m_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x[m_col - 1] = b.m_elem[m_row - 1];
    for (k = m_col - 2, pU += m_row * (m_col - 2); k >= 0; pU -= m_col, k--)
    {
        x[k] = b.m_elem[k];
        for (i = k + 1; i < m_col; i++)
            x[k] -= x[i] * *(pU + i);
    }

    if (isOk)
        *isOk = 1;

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtUpperTriangular<m_row, m_col, m_type>::InverseUnit(dtMatrix<m_row, m_col, m_type> U, dtMatrix<m_row, m_col, m_type> &invU)
{
    if (m_row != m_col)
        return -1;

    int i, j, k;
    m_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = m_row - 2, pMi = U.m_elem + m_col * (m_row - 2); i >= 0; pMi -= m_col, i--)
    {
        for (j = m_col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + m_col; k < j; pMk += m_col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invU = U;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtUpperTriangular<m_row, m_col, m_type>::InverseUnit(dtMatrix<m_row, m_col, m_type> U, int8_t *isOk)
{
    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtMatrix<m_row, m_col, m_type>();
    }

    int i, j, k;
    m_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = m_row - 2, pMi = U.m_elem + m_col * (m_row - 2); i >= 0; pMi -= m_col, i--)
    {
        for (j = m_col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + m_col; k < j; pMk += m_col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    if (isOk)
        *isOk = 1;

    return U;
}

} // namespace dtMath

#endif // DTMATH_DTUPPER_TRIANGULAR_TPP_