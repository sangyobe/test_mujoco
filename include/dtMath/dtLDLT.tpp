/*!
\file       dtLDLT.h
\brief      dtMath, Cholesky decomposition(L*D*L^T form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLDLT_TPP_
#define DTMATH_DTLDLT_TPP_

#include "dtLDLT.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLDLT<m_row, m_col, m_type>::dtLDLT()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLDLT<m_row, m_col, m_type>::dtLDLT(const m_type *element, const size_t n_byte)
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
inline dtLDLT<m_row, m_col, m_type>::dtLDLT(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLDLT<m_row, m_col, m_type>::dtLDLT(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Compute()
{
    if (m_row != m_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type Lik;

    for (i = 1, pMi = m_elem + m_col; i < m_row; pMi += m_col, i++)
    {
        // Calculate the part of L*D Matrix from Lij*Djj
        // L*D matrix are used for Calculating Dii
        // here *(pMi + j) = Lij*Djj
        for (j = 0, pMj = m_elem; j < i; j++, pMj += m_row)
            for (k = 0; k < j; k++)
                *(pMi + j) -= *(pMi + k) * *(pMj + k);

        // Calculate the diagonal element Dii
        // Calculate the Mij from Lik and Store the Mjk
        for (k = 0, pMk = m_elem; k < i; pMk += m_col, k++)
        {
            Lik = *(pMi + k) / *(pMk + k);  // *(pMi + k) = Lik * Dkk
            *(pMi + i) -= *(pMi + k) * Lik; // Dii = Aii - sum(Lik * Dkk * Lik)
            *(pMi + k) = Lik;               // Lik
            *(pMk + i) = Lik;               // Lki
        }

        // If diagonal element is not positive, return the error,
        // the matrix is not positive definite symmetric.
        if (*(pMi + i) <= std::numeric_limits<m_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
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
inline int8_t dtLDLT<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLDLT<m_row, m_col, m_type>::GetMatrix() const
{
    return dtMatrix<m_row, m_col, m_type>(m_elem);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLDLT<m_row, m_col, m_type>::GetMatrixL() const
{
    int i, j;
    m_type L[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        L[i * m_col + i] = 1;

    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            L[i * m_col + j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(L);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLDLT<m_row, m_col, m_type>::GetMatrixD() const
{
    int i;
    m_type D[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        D[i * m_col + i] = m_elem[i * m_col + i];

    return dtMatrix<m_row, m_col, m_type>(D);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLDLT<m_row, m_col, m_type>::GetMatrixU() const
{
    int i, j;
    m_type U[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        U[i * m_col + i] = 1;

    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            U[i + m_col * j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(U);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline int8_t dtLDLT<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, dtMatrix<m_col, col, m_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    m_type *pMi;

    if (!m_isOk)
        return -1;

    for (j = 0; j < col; j++)
    {
        /* Solve Lz = b */
        // Solve the unit lower triangular (forward substitution), here x is z
        for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
        {
            x.m_elem[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];
        }

        /* Solve Dy = z */
        // Solve the diagonal matrix, here x is y
        for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
        {
            x.m_elem[i * col + j] /= *(pMi + i);
        }

        /* Solve Ux = y */
        // Solve the unit upper triangular (backward substitution)
        for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
        {
            for (k = i + 1; k < m_col; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;

    if (!m_isOk)
        return -1;

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline dtMatrix<m_col, col, m_type> dtLDLT<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    m_type *pMi;
    m_type x[m_col * col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtMatrix<m_col, col, m_type>(x);
    }

    for (j = 0; j < col; j++)
    {
        /* Solve Lz = b */
        // Solve the unit lower triangular (forward substitution), here x is z
        for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
        {
            x[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];
        }

        /* Solve Dy = z */
        // Solve the diagonal matrix (divided by diagonal elements), here x is y
        for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
        {
            x[i * col + j] /= *(pMi + i);
        }

        /* Solve Ux = y */
        // Solve the upper triangular (backward substitution)
        for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
        {
            for (k = i + 1; k < m_col; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];
        }
    }

    return dtMatrix<m_col, col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtLDLT<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
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
        return dtVector<m_col, m_type>(x);
    }

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix (divided by diagonal elements), here x is y
    for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
    {
        x[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> &inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::Inverse(dtMatrix3<m_type, m_row, m_col> &inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLDLT<m_row, m_col, m_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type inv[m_row * m_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return dtMatrix<m_row, m_col, m_type>();
    }

    memcpy(inv, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return dtMatrix<m_row, m_col, m_type>(inv);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLDLT<m_row, m_col, m_type>::InverseArray(m_type *inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type *dtLDLT<m_row, m_col, m_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type inv[m_row * m_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return inv;
    }

    memcpy(inv, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return inv;
}

} // namespace dtMath

#endif // DTMATH_DTLDLT_TPP_
