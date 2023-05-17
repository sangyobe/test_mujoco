/*!
\file       dtHouseholderQR.h
\brief      dtMath, QR Decomposition without pivoting(Householder method) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQR_TPP_
#define DTMATH_DTQR_TPP_

#include "dtQR.h"

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtQR<m_row, m_col, m_type>::dtQR()
{
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    m_isOk = 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtQR<m_row, m_col, m_type>::dtQR(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
        m_isOk = 0;
    else
    {
        memcpy(m_elem, element, n_byte);
        Compute();
    }
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtQR<m_row, m_col, m_type>::dtQR(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtQR<m_row, m_col, m_type>::dtQR(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtQR<m_row, m_col, m_type>::Compute()
{
    // A = Q*R
    int col = (m_row > m_col) ? m_col : m_row;
    int i, j, k, jj;

    m_type H1n[m_row * m_row]; // size: mxm
    m_type Hn1[m_row * m_col] = {
        0,
    }; // Hn*...*H1*A, mxn
    m_type H[m_row * m_row];
    m_type u[m_row];
    m_type *pQi, *pQij;
    m_type *pRi, *pRij;
    m_type *pHi, *pHkj;
    m_type *pH1ni;
    m_type *pHn1kj;
    m_type *pUi, *pUij;
    m_type sqNorm;

    /* QR Decomposition */
    // R = Hn * ... * H1 * A, size:mxn
    // Q = H1 * ... * Hn, size:mxm
    // H = I - 2 * u * uT / uT * u, size:mxm
    // u = a + sign(a1) * norm(a) * e1
    // a is column vector of R or A
    // e1 = [1 0 .. 0]T

    // copy mat A to mat Hn1A
    memcpy(Hn1, m_elem, sizeof(m_type) * m_row * m_col);

    for (j = 0, pUi = Hn1; j < col; j++, pUi += m_col)
    {
        // calculate the vector u
        sqNorm = 0;
        for (i = j, pUij = pUi + j; i < m_row; i++, pUij += m_col)
        {
            u[i] = *(pUij);              // vector a (column of R or Hn1*A)
            sqNorm += *(pUij) * *(pUij); // squared norm of vector a
        }

        u[j] += (u[j] > 0) ? std::sqrt(sqNorm) : -std::sqrt(sqNorm);

        // uT*u, squared norm of vector u
        sqNorm = 0;
        for (i = j; i < m_row; i++)
            sqNorm += *(u + i) * *(u + i);

        if (sqNorm > std::numeric_limits<m_type>::epsilon())
        {
            // calculate the householder matrix H = I - 2*u*uT / uT*u
            for (i = j, pHi = H + j * m_row; i < m_row; i++, pHi += m_row)
            {
                for (k = j; k < m_row; k++)
                {
                    if (i == k)
                        *(pHi + k) = 1 - 2 * *(u + i) * *(u + k) / sqNorm;
                    else
                        *(pHi + k) = -2 * *(u + i) * *(u + k) / sqNorm;
                }
            }

            if (j == 0)
            {
                // Q = I*H1;
                memcpy(H1n, H, sizeof(m_type) * m_row * m_row);
                memcpy(m_Q, H, sizeof(m_type) * m_row * m_row);

                // R = H1*A(=Hn1)
                pRi = m_R;
                pHi = H;
                for (i = 0; i < m_row; i++, pRi += m_col, pHi += m_row)
                {
                    for (jj = 0; jj < m_col; jj++)
                    {
                        pHn1kj = Hn1 + jj;
                        pRij = pRi + jj;
                        for (k = 0; k < m_row; k++, pHn1kj += m_col)
                            *(pRij) += *(pHi + k) * *pHn1kj;
                    }
                }
                memcpy(Hn1, m_R, sizeof(m_type) * m_row * m_col);
            }
            else
            {
                // Q = (H1*...*Hn-1) * H
                pQi = m_Q;
                pH1ni = H1n;
                for (i = 0; i < m_row; i++, pQi += m_row, pH1ni += m_row)
                {
                    for (jj = j; jj < m_row; jj++)
                    {
                        pHkj = H + m_row * j + jj;
                        pQij = pQi + jj;
                        *pQij = 0;
                        for (k = j; k < m_row; k++, pHkj += m_row)
                            *(pQij) += *(pH1ni + k) * *pHkj;
                    }
                }
                memcpy(H1n, m_Q, sizeof(m_type) * m_row * m_row);

                // R = H * (Hn-1*...*H1*A)
                pRi = m_R + m_col * j;
                pHi = H + m_row * j;
                for (i = j; i < m_row; i++, pRi += m_col, pHi += m_row)
                {
                    for (jj = j; jj < m_col; jj++)
                    {
                        pHn1kj = Hn1 + m_col * j + jj;
                        pRij = pRi + jj;
                        *pRij = 0;
                        for (k = j; k < m_row; k++, pHn1kj += m_col)
                            *pRij += *(pHi + k) * *pHn1kj;
                    }
                }
                memcpy(Hn1, m_R, sizeof(m_type) * m_row * m_col);
            }
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtQR<m_row, m_col, m_type>::Compute(const m_type *element, const size_t n_byte)
{
    if ((sizeof(m_type) * m_row * m_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memcpy(m_elem, element, n_byte);
    memset(m_Q, 0, sizeof(m_type) * m_row * m_row);
    memset(m_R, 0, sizeof(m_type) * m_row * m_col);

    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtQR<m_row, m_col, m_type>::Compute(const dtMatrix<m_row, m_col, m_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    memset(m_Q, 0, sizeof(m_type) * m_row * m_row);
    memset(m_R, 0, sizeof(m_type) * m_row * m_col);

    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t dtQR<m_row, m_col, m_type>::Compute(const dtMatrix3<m_type, m_row, m_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(m_type) * m_row * m_col);
    memset(m_Q, 0, sizeof(m_type) * m_row * m_row);
    memset(m_R, 0, sizeof(m_type) * m_row * m_col);

    return Compute();
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_row, m_type> dtQR<m_row, m_col, m_type>::GetMatrixQ() const
{
    return dtMatrix<m_row, m_row, m_type>(m_Q);
}

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline dtMatrix<m_row, m_col, m_type> dtQR<m_row, m_col, m_type>::GetMatrixR() const
{
    return dtMatrix<m_row, m_col, m_type>(m_R);
}

} // namespace dtMath

#endif // DTMATH_DTQR_TPP_