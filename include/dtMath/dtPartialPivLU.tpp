/*!
\file       dtPartialPivLU.h
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTPARTIAL_PIV_LU_TPP_
#define DTMATH_DTPARTIAL_PIV_LU_TPP_

#include "dtPartialPivLU.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtPartialPivLU<m_row, m_col, m_type>::dtPartialPivLU()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtPartialPivLU<m_row, m_col, m_type>::dtPartialPivLU(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
        m_isOk = 0;
    else
    {
        memset(m_elem, element, n_byte);
        Compute();
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtPartialPivLU<m_row, m_col, m_type>::dtPartialPivLU(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtPartialPivLU<m_row, m_col, m_type>::dtPartialPivLU(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Compute()
{
    if (m_row != m_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    m_type *pMi, *pMk, *p_pivotRow = nullptr;
    m_type max, absElem;
    m_type pivotRow[m_col];

    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        /* Pivoting */
        // find the pivot row
        m_pivot[i] = i;
        max = std::abs(*(pMi + i));
        for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
        {
            if (max < (absElem = std::abs(*(pMk + i))))
            {
                max = absElem;
                m_pivot[i] = k;
                p_pivotRow = pMk; // pMk is pivot row
            }
        }

        // interchange the two rows.
        if (m_pivot[i] != i)
        {
            memcpy(pivotRow, p_pivotRow, sizeof(m_type) * m_col);
            memcpy(p_pivotRow, pMi, sizeof(m_type) * m_col);
            memcpy(pMi, pivotRow, sizeof(m_type) * m_col);
        }

        // matrix is singular, return error
        if (std::abs(*(pMi + i)) <= std::numeric_limits<m_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        /* LU Decompostion using Gaussian Elimination */
        for (k = i + 1, pMk = pMi + m_col; k < m_row; pMk += m_col, k++)
        {
            // find the lower triangular matrix elements for column i.
            *(pMk + i) /= *(pMi + i);

            // update the upper triangular matrix for remaining matrix
            for (j = i + 1; j < m_col; j++)
                *(pMk + j) -= *(pMk + i) * *(pMi + j);
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memset(m_elem, element, n_byte);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type dtPartialPivLU<m_row, m_col, m_type>::Determinant()
{
    if (!m_isOk)
        return -1;
    uint16_t offset = m_row + 1;
    m_type det = 1;

    for (uint16_t i = 0; i < m_row; i++)
        det *= m_elem[i * offset];

    return det;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::GetMatrix() const
{
    return dtMatrix<m_row, m_col, m_type>(m_elem);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::GetMatrixL() const
{
    int i, j;
    m_type L[m_row * m_col] = {
        0,
    };

    /* Set diagonal elements as 1 */
    for (i = 0; i < m_row; i++)
    {
        L[i * (m_col + 1)] = 1;
    }

    /* Update remaining matrix from m_elem to L*/
    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            L[i * m_col + j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(L);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::GetMatrixU() const
{
    int i, j;
    m_type U[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        for (j = i; j < m_col; j++)
            U[i * m_col + j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(U);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::GetMatrixP() const
{
    int i;
    m_type P[m_row * m_col] = {
        0,
    };
    m_type row[m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        P[i * (m_col + 1)] = 1;

    for (i = 0; i < m_row; i++)
    {
        if (m_pivot[i] != i)
        {
            memcpy(row, &P[i * m_col], sizeof(m_type) * m_col);
            memcpy(&P[i * m_col], &P[m_pivot[i] * m_col], sizeof(m_type) * m_col);
            memcpy(&P[m_pivot[i] * m_col], row, sizeof(m_type) * m_col);
        }
    }

    return dtMatrix<m_row, m_col, m_type>(P);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;
    m_type tmp;
    m_type vb[m_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    if (!m_isOk)
        return -1;

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x.m_elem[i] = vb[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<m_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;
    m_type x[m_col] = {
        0,
    };
    m_type tmp;
    m_type vb[m_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtVector<m_col, m_type>();
    }

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x[i] = vb[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < m_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<m_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> &inv)
{
    int i, j, k;
    int colIdx[m_col]; // column index for interchange the current col with the pivot col
    m_type *p_Mi;
    m_type *p_invMi, *p_invMj, *p_invMk;
    m_type sum;
    m_type invL[m_row * m_col] = {
        0,
    };
    m_type invU[m_row * m_col] = {
        0,
    };

    if (!m_isOk)
        return -1;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                colIdx[j] = colIdx[m_pivot[j]];
                colIdx[m_pivot[j]] = j;
            }

            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::Inverse(dtMatrix3<m_type, m_row, m_col> &inv)
{
    int i, j, k;
    int colIdx[m_col]; // column index for interchange the current col with the pivot col
    m_type *p_Mi;
    m_type *p_invMi, *p_invMj, *p_invMk;
    m_type sum;
    m_type invL[m_row * m_col] = {
        0,
    };
    m_type invU[m_row * m_col] = {
        0,
    };

    if (!m_isOk)
        return -1;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                colIdx[j] = colIdx[m_pivot[j]];
                colIdx[m_pivot[j]] = j;
            }

            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtPartialPivLU<m_row, m_col, m_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    int colIdx[m_col]; // column index for interchange the current col with the pivot col
    m_type *p_Mi;
    m_type *p_invMi, *p_invMj, *p_invMk;
    m_type sum;
    m_type invL[m_row * m_col] = {
        0,
    };
    m_type invU[m_row * m_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtMatrix<m_row, m_col, m_type>();
    }

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    memset(m_inv, 0, sizeof(m_type) * m_row * m_col);
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                colIdx[j] = colIdx[m_pivot[j]];
                colIdx[m_pivot[j]] = j;
            }

            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return dtMatrix<m_row, m_col, m_type>(m_inv);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtPartialPivLU<m_row, m_col, m_type>::InverseArray(m_type *inv)
{
    int i, j, k;
    int colIdx[m_col]; // column index for interchange the current col with the pivot col
    m_type *p_Mi;
    m_type *p_invMi, *p_invMj, *p_invMk;
    m_type sum;
    m_type invL[m_row * m_col] = {
        0,
    };
    m_type invU[m_row * m_col] = {
        0,
    };

    if (!m_isOk)
        return -1;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                colIdx[j] = colIdx[m_pivot[j]];
                colIdx[m_pivot[j]] = j;
            }

            for (k = 0; k < m_col; k++)
                inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type *dtPartialPivLU<m_row, m_col, m_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    int colIdx[m_col]; // column index for interchange the current col with the pivot col
    m_type *p_Mi;
    m_type *p_invMi, *p_invMj, *p_invMk;
    m_type sum;
    m_type invL[m_row * m_col] = {
        0,
    };
    m_type invU[m_row * m_col] = {
        0,
    };
    memset(m_inv, 0, sizeof(m_type) * m_row * m_col);

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return m_inv;
    }

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                colIdx[j] = colIdx[m_pivot[j]];
                colIdx[m_pivot[j]] = j;
            }

            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return m_inv;
}

} // namespace dtMath

#endif // DTMATH_DTPARTIAL_PIV_LU_TPP_