/*!
\file       dtLowerTriangular.h
\brief      dtMath, Lower triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLOWER_TRIANGULAR_TPP_
#define DTMATH_DTLOWER_TRIANGULAR_TPP_

#include "dtLowerTriangular.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLowerTriangular<m_row, m_col, m_type>::Solve(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    if (m_row != m_col)
        return -1;

    m_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < m_row; i++, pMi += m_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<m_type>::epsilon())
            return -1; // singular;

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtLowerTriangular<m_row, m_col, m_type>::Solve(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, int8_t *isOk)
{
    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    m_type x[m_col] = {
        0,
    };
    m_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < m_row; i++, pMi += m_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<m_type>::epsilon())
        {
            if (isOk)
                *isOk = 0; // singular;
            return dtVector<m_col, m_type>();
        }

        x[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    if (isOk)
        *isOk = 1;

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLowerTriangular<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> L, dtMatrix<m_row, m_col, m_type> &invL)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;

    if (m_row != m_col)
        return -1;

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<m_type>::epsilon())
            return -1;
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invL = L;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLowerTriangular<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> L, int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;

    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtMatrix<m_row, m_col, m_type>();
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        if (*pMk == 0) // To do dhl
        {
            if (isOk)
                *isOk = 0;
            return dtMatrix<m_row, m_col, m_type>();
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    if (isOk)
        *isOk = 1;

    return L;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLowerTriangular<m_row, m_col, m_type>::SolveUnit(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    if (m_row != m_col)
        return -1;

    int i, k;
    m_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += m_col; k < m_row; pL += m_col, k++)
    {
        x.m_elem[i] = b.m_elem[i];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtLowerTriangular<m_row, m_col, m_type>::SolveUnit(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, int8_t *isOk)
{
    if (m_row != m_col)
    {
        if (isOk)
            *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    int i, k;
    m_type x[m_col] = {
        0,
    };
    m_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x[0] = b.m_elem[0];
    for (k = 1, pL += m_col; k < m_row; pL += m_col, k++)
    {
        x[i] = b.m_elem[i];
        for (i = 0, x[k] = b.m_elem[k]; i < k; i++)
            x[k] -= x[i] * *(pL + i);
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLowerTriangular<m_row, m_col, m_type>::InverseUnit(dtMatrix<m_row, m_col, m_type> L, dtMatrix<m_row, m_col, m_type> &invL)
{
    if (m_row != m_col)
        return -1;

    int i, j, k;
    m_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invL = L;

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLowerTriangular<m_row, m_col, m_type>::InverseUnit(dtMatrix<m_row, m_col, m_type> L, int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    if (isOk)
        *isOk = 1;

    return L;
}

} // namespace dtMath

#endif // DTMATH_DTLOWER_TRIANGULAR_TPP_
