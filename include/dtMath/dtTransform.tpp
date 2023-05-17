/*!
\file       dtTransform.h
\brief      dtMath, Homogeneous transformation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTTRANSFORM_TPP_
#define DTMATH_DTTRANSFORM_TPP_

#include "dtTransform.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col>::dtTransform()
{
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col>::dtTransform(const dtRotation<m_type, 3, 3> &R, const dtVector3<m_type, 3> &p)
{
    m_R = R;
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col>::dtTransform(const dtQuaternion<m_type, 4> &q, const dtVector3<m_type, 3> &p)
{
    m_R = dtRotation<m_type, 3, 3>(q);
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col>::dtTransform(const uint16_t order, const dtVector3<m_type, 3> &e, const dtVector3<m_type, 3> &p)
{
    m_R = dtRotation<m_type, 3, 3>(order, e);
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col>::dtTransform(const dtTransform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetZero()
{
    m_R.SetZero();
    m_p.SetZero();
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetIdentity()
{
    m_R.SetIdentity();
    m_p.SetZero();
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtVector3<m_type, 3> &p)
{
    m_p = p;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtRotation<m_type, 3, 3> &R)
{
    m_R = R;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtQuaternion<m_type, 4> &q)
{
    m_R = dtRotation<m_type, 3, 3>(q);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const uint16_t order, const dtVector3<m_type, 3> &e)
{
    m_R = dtRotation<m_type, 3, 3>(order, e);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtRotation<m_type, 3, 3> &R, const dtVector3<m_type, 3> &p)
{
    m_R = R;
    m_p = p;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtQuaternion<m_type, 4> &q, const dtVector3<m_type, 3> &p)
{
    m_R = dtRotation<m_type, 3, 3>(q);
    m_p = p;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const uint16_t order, const dtVector3<m_type, 3> &e, const dtVector3<m_type, 3> &p)
{
    m_R = dtRotation<m_type, 3, 3>(order, e);
    m_p = p;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::SetElement(const dtTransform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector6<m_type, 6> dtTransform<m_type, m_row, m_col>::GetError(const dtTransform &m) const
{
    // dtVector3<m_type, 3> pos_error = m_p - m.m_p;
    // dtQuaternion<m_type, 4> qd(m_R);
    // dtVector3<m_type, 3> ori_error = qd.GetOriErr(dtQuaternion<m_type, 4>(m.m_R));
    //
    // return dtVector6<m_type, 6>(pos_error, ori_error);

    return dtVector6<m_type, 6>(
        m_p - m.m_p,
        dtQuaternion<m_type, 4>(m_R).GetOriErr(dtQuaternion<m_type, 4>(m.m_R)));
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix<m_col, m_row, m_type> dtTransform<m_type, m_row, m_col>::Transpose() const
{
    m_type mat[m_col * m_row] = {
        m_R.m_elem[0], m_R.m_elem[3], m_R.m_elem[6], 0,
        m_R.m_elem[1], m_R.m_elem[4], m_R.m_elem[7], 0,
        m_R.m_elem[2], m_R.m_elem[5], m_R.m_elem[8], 0,
        m_p.m_elem[0], m_p.m_elem[1], m_p.m_elem[2], (m_type)(1)};

    return dtMatrix<m_col, m_row, m_type>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col> dtTransform<m_type, m_row, m_col>::Inv() const
{
    return dtTransform(m_R.Inv(), -(m_R.Inv() * m_p));
}

/* Member access operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline m_type &dtTransform<m_type, m_row, m_col>::operator()(uint16_t irow, uint16_t icol)
{
    if (irow < 3 && icol == 3)
        return m_p.m_elem[irow];
    else if (irow == 3)
        return m_dummy[icol];

    return m_R.m_elem[irow * 3 + icol];
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline const m_type &dtTransform<m_type, m_row, m_col>::operator()(uint16_t irow, uint16_t icol) const
{
    if (irow < 3 && icol == 3)
        return m_p.m_elem[irow];
    else if (irow == 3)
        return m_dummy[icol];

    return m_R.m_elem[irow * 3 + icol];
}

/* Assignment operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col> &dtTransform<m_type, m_row, m_col>::operator=(const dtTransform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;

    return (*this);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix<m_row, m_col, m_type> dtTransform<m_type, m_row, m_col>::operator+(const dtMatrix<m_row, m_col, m_type> &m) const
{
    m_type mat[m_row * m_col] = {
        m_R.m_elem[0] + m.m_elem[0], m_R.m_elem[1] + m.m_elem[1], m_R.m_elem[2] + m.m_elem[2], m_p.m_elem[0] + m.m_elem[3],
        m_R.m_elem[3] + m.m_elem[4], m_R.m_elem[4] + m.m_elem[5], m_R.m_elem[5] + m.m_elem[6], m_p.m_elem[1] + m.m_elem[7],
        m_R.m_elem[6] + m.m_elem[8], m_R.m_elem[7] + m.m_elem[9], m_R.m_elem[8] + m.m_elem[10], m_p.m_elem[2] + m.m_elem[11],
        m.m_elem[12], m.m_elem[13], m.m_elem[14], (m_type)(1) + m.m_elem[15]};

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtMatrix<m_row, m_col, m_type> dtTransform<m_type, m_row, m_col>::operator-(const dtMatrix<m_row, m_col, m_type> &m) const
{
    m_type mat[m_row * m_col] = {
        m_R.m_elem[0] - m.m_elem[0], m_R.m_elem[1] - m.m_elem[1], m_R.m_elem[2] - m.m_elem[2], m_p.m_elem[0] - m.m_elem[3],
        m_R.m_elem[3] - m.m_elem[4], m_R.m_elem[4] - m.m_elem[5], m_R.m_elem[5] - m.m_elem[6], m_p.m_elem[1] - m.m_elem[7],
        m_R.m_elem[6] - m.m_elem[8], m_R.m_elem[7] - m.m_elem[9], m_R.m_elem[8] - m.m_elem[10], m_p.m_elem[2] - m.m_elem[11],
        -m.m_elem[12], -m.m_elem[13], -m.m_elem[14], (m_type)(1) - m.m_elem[15]};

    return dtMatrix<m_row, m_col, m_type>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
template <uint16_t col>
inline dtMatrix<m_row, col, m_type> dtTransform<m_type, m_row, m_col>::operator*(const dtMatrix<m_row, col, m_type> &m) const
{
    m_type mat[m_row * col];

    for (uint16_t irow = 0; irow < 3; ++irow)
    {
        for (uint16_t icol = 0; icol < col; ++icol)
        {
            mat[irow * col + icol] =
                m_R.m_elem[3 * irow] * m.m_elem[icol] +
                m_R.m_elem[3 * irow + 1] * m.m_elem[col + icol] +
                m_R.m_elem[3 * irow + 2] * m.m_elem[2 * col + icol] +
                m_p.m_elem[irow] * m.m_elem[3 * col + icol];
        }
    }

    memcpy(&mat[3 * col], &m.m_elem[3 * col], sizeof(m_type) * col);

    return dtMatrix<m_row, col, m_type>(mat);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtTransform<m_type, m_row, m_col> dtTransform<m_type, m_row, m_col>::operator*(const dtTransform &m) const
{
    return dtTransform(m_R * m.m_R, m_R * m.m_p + m_p);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtTransform<m_type, m_row, m_col>::operator*(const dtVector<3, m_type> &v) const
{
    return dtVector3<m_type, 3>(m_R * v + m_p);
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline dtVector3<m_type, 3> dtTransform<m_type, m_row, m_col>::operator*(const dtVector3<m_type, 3> &v) const
{
    return dtVector3<m_type, 3>(m_R * v + m_p);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtTransform<m_type, m_row, m_col>::operator==(const dtTransform &m) const
{
    if (m_R != m.m_R)
        return false;
    if (m_p != m.m_p)
        return false;

    return true;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline bool dtTransform<m_type, m_row, m_col>::operator!=(const dtTransform &m) const
{
    if (m_R != m.m_R)
        return true;
    if (m_p != m.m_p)
        return true;

    return false;
}

template <typename m_type, uint16_t m_row, uint16_t m_col>
inline void dtTransform<m_type, m_row, m_col>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < 3; irow++)
    {
        Serial.printf("%7.3f %7.3f %7.3f %7.3f\n",
                      m_R.m_elem[irow * 3], m_R.m_elem[irow * 3 + 1], m_R.m_elem[irow * 3 + 2], m_p.m_elem[irow]);
    }
    Serial.printf("%7.3f %7.3f %7.3f %7.3f\n", m_dummy[0], m_dummy[1], m_dummy[2], m_dummy[3]);
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < 3; irow++)
    {
        printf("%7.3f %7.3f %7.3f %7.3f\n",
               m_R.m_elem[irow * 3], m_R.m_elem[irow * 3 + 1], m_R.m_elem[irow * 3 + 2], m_p.m_elem[irow]);
    }
    printf("%7.3f %7.3f %7.3f %7.3f\n", m_dummy[0], m_dummy[1], m_dummy[2], m_dummy[3]);
    printf("%c", endChar);
#endif
}

typedef dtTransform<> dtTMat;

} // namespace dtMath

#endif // DTMATH_DTTRANSFORM_TPP_
