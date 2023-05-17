/*!
\file       dtNoPivLU.h
\brief      dtMath, LU Decomposition without pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTNO_PIV_LU_TPP_
#define DTMATH_DTNO_PIV_LU_TPP_

#include "dtNoPivLU.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtNoPivLU<m_row, m_col, m_type>::dtNoPivLU()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtNoPivLU<m_row, m_col, m_type>::dtNoPivLU(const m_type *element, const size_t n_byte)
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
inline dtNoPivLU<m_row, m_col, m_type>::dtNoPivLU(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtNoPivLU<m_row, m_col, m_type>::dtNoPivLU(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Compute()
{
    if (m_row != m_col)
    {
        m_isOk = 0;
        return -1;
    }

    int x, i, j, k;
    m_type *pMx, *pMi, *pMk;

    for (x = 0, pMx = m_elem; x < m_row; pMx += m_col, x++)
    {
        /* find the uppper triangular matrix element for row x */
        for (j = x; j < m_col; j++)
        {
            // U(x,j) = A(x,j) - sum of (L(x,k) * U(k,j))
            for (k = 0, pMk = m_elem; k < x; pMk += m_col, k++)
                *(pMx + j) -= *(pMx + k) * *(pMk + j); // in-place
        }

        if (std::abs(*(pMx + x)) <= std::numeric_limits<m_type>::epsilon())
        {
            m_isOk = 0;
            return -1; // matrix is singular, if Diagonal is 0
        }

        /* find the lower triangular matrix element for col x */
        for (i = x + 1, pMi = pMx + m_row; i < m_row; pMi += m_row, i++)
        {
            // L(i,x) = [A(i,x) - sum of (L(i,k) * U(k,x))] / Uxx
            for (k = 0, pMk = m_elem; k < x; pMk += m_row, k++)
                *(pMi + x) -= *(pMi + k) * *(pMk + x);
            *(pMi + x) /= *(pMx + x);
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
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
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtNoPivLU<m_row, m_col, m_type>::GetMatrix() const
{
    return dtMatrix<m_row, m_col, m_type>(m_elem);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtNoPivLU<m_row, m_col, m_type>::GetMatrixL() const
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
inline dtMatrix<m_row, m_col, m_type> dtNoPivLU<m_row, m_col, m_type>::GetMatrixU() const
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
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;

    if (!m_isOk)
        return -1;

    /* Solve Ly = b */
    // Solve the unit lower triangular (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + k)) <= std::numeric_limits<m_type>::epsilon()) return -1; // The matrix U is singular

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtNoPivLU<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;
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

    /* Solve Ly = b */
    // Unit Lower Triangular Solve (forward substitution)
    for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
    {
        x[i] = b.m_elem[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // Upper Triangular Solve (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pU + k)) <= std::numeric_limits<m_type>::epsilon()) m_isOk = 0; // The matrix U is singular

        x[i] /= *(pMi + i);
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> &inv)
{
    int i, j, k;
    m_type *pMi;
    m_type *pInvMi, *pInvMj, *pInvMk;
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
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::Inverse(dtMatrix3<m_type, m_row, m_col> &inv)
{
    int i, j, k;
    m_type *pMi;
    m_type *pInvMi, *pInvMj, *pInvMk;
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
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtNoPivLU<m_row, m_col, m_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi;
    m_type *pInvMi, *pInvMj, *pInvMk;
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
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    memset(m_inv, 0, sizeof(m_type) * m_row * m_col);
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(m_inv);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtNoPivLU<m_row, m_col, m_type>::InverseArray(m_type *inv)
{
    int i, j, k;
    m_type *pMi;
    m_type *pInvMi, *pInvMj, *pInvMk;
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
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                inv[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type *dtNoPivLU<m_row, m_col, m_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi;
    m_type *pInvMi, *pInvMj, *pInvMk;
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
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<m_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return m_inv;
}

} // namespace dtMath

#endif // DTMATH_DTNO_PIV_LU_TPP_