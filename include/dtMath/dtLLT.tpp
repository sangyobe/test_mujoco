/*!
\file       dtLLT.h
\brief      dtMath, Cholesky decomposition(L*L^T form) Class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLLT_TPP_
#define DTMATH_DTLLT_TPP_

#include "dtLLT.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLLT<m_row, m_col, m_type>::dtLLT()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLLT<m_row, m_col, m_type>::dtLLT(const m_type *element, const size_t n_byte)
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
inline dtLLT<m_row, m_col, m_type>::dtLLT(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtLLT<m_row, m_col, m_type>::dtLLT(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Compute()
{
    if (m_row != m_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    m_type *pMi, *pMj, *pMjj, *pMjk;

    for (j = 0, pMj = m_elem; j < m_col; pMj += m_row, j++)
    {
        /* Calculate the Diagonal element in colum j */
        pMjj = pMj + j;
        for (k = 0, pMjk = pMj; k < j; pMjk += 1, k++)
            *pMjj -= *pMjk * *pMjk;

        // If diagonal element is not positive, return the error,
        // the matrix is not positive definite symmetric.
        if (*pMjj <= std::numeric_limits<m_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        *pMjj = std::sqrt(*pMjj);

        /* Calculate the lower triangular matrix for colum j */
        for (i = j + 1, pMi = pMj + m_col; i < m_row; pMi += m_col, i++)
        {
            for (k = 0; k < j; k++)
                *(pMi + j) -= *(pMi + k) * *(pMj + k);

            *(pMi + j) /= *pMjj;     // Lower Triangular Matrix, in-place
            *(pMj + i) = *(pMi + j); // Upper Triangular Matrix, in-place
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
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
inline int8_t dtLLT<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLLT<m_row, m_col, m_type>::GetMatrix() const
{
    return dtMatrix<m_row, m_col, m_type>(m_elem);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLLT<m_row, m_col, m_type>::GetMatrixL() const
{
    int i, j;
    m_type L[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        for (j = 0; j <= i; j++)
            L[i * m_col + j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(L);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLLT<m_row, m_col, m_type>::GetMatrixU() const
{
    int i, j;
    m_type U[m_row * m_col] = {
        0,
    };

    for (i = 0; i < m_row; i++)
        for (j = 0; j < i; j++)
            U[i + m_col * j] = m_elem[i * m_col + j];

    return dtMatrix<m_row, m_col, m_type>(U);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline int8_t dtLLT<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, dtMatrix<m_col, col, m_type> &x)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    m_type *pMi;

    if (!m_isOk)
        return -1;

    for (j = 0; j < col; j++)
    {
        /* Solve Ly = b */
        // Solve the lower triangular matrix for y (forward substitution), here x is y
        for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
        {
            x.m_elem[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];

            x.m_elem[i * col + j] /= *(pMi + i);
        }

        /* Solve LTx = y */
        // Solve the upper triangular (backward substitution), LT = U
        for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
        {
            for (k = i + 1; k < m_col; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];

            x.m_elem[i * col + j] /= *(pMi + i);
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    m_type *pMi;

    if (!m_isOk)
        return -1;

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve LTx = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
template <uint16_t col>
inline dtMatrix<m_col, col, m_type> dtLLT<m_row, m_col, m_type>::Solve(const dtMatrix<m_row, col, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
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
        /* Solve Ly = b */
        // Solve the lower triangular matrix for y (forward substitution), here x is y
        for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
        {
            x[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];

            x[i * col + j] /= *(pMi + i);
        }

        /* Solve LTx = y */
        // Solve the upper triangular (backward substitution), LT = U
        for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
        {
            for (k = i + 1; k < m_col; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];

            x[i * col + j] /= *(pMi + i);
        }
    }

    return dtMatrix<m_col, col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtVector<m_col, m_type> dtLLT<m_row, m_col, m_type>::Solve(const dtVector<m_row, m_type> &b, int8_t *isOk)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
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
        return dtVector<m_col, m_type>(x);
    }

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    /* Solve LTx = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    return dtVector<m_col, m_type>(x);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Inverse(dtMatrix<m_row, m_col, m_type> &inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv.m_elem; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j <= i; j++, pMj += m_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < m_col; k++, pMk += m_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::Inverse(dtMatrix3<m_type, m_row, m_col> &inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv.m_elem; k < m_row; pMk += (m_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv.m_elem; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j <= i; j++, pMj += m_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < m_col; k++, pMk += m_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtLLT<m_row, m_col, m_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;
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

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < m_row; pMk += (m_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += m_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < m_col; k++, pMk += m_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return dtMatrix<m_row, m_col, m_type>(inv);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtLLT<m_row, m_col, m_type>::InverseArray(m_type *inv)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv, m_elem, sizeof(m_type) * m_row * m_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < m_row; pMk += (m_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += m_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < m_col; k++, pMk += m_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline m_type *dtLLT<m_row, m_col, m_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    m_type *pMi, *pMj, *pMk;
    m_type sum;
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

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < m_row; pMk += (m_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += m_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += m_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < m_col; k++, pMk += m_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return inv;
}

} // namespace dtMath

#endif // DTMATH_DTLLT_TPP_
