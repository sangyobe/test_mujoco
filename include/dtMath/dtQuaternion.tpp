/*!
\file       dtQuaternion.h
\brief      dtMath, Quaternion class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQUATERNION_TPP_
#define DTMATH_DTQUATERNION_TPP_

#include "dtQuaternion.h"

namespace dtMath
{

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(/* args */)
{
    m_elem[0] = 1;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const m_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const m_type w, const m_type x, const m_type y, const m_type z)
{
    m_elem[0] = w;
    m_elem[1] = x;
    m_elem[2] = y;
    m_elem[3] = z;
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const uint16_t order, const m_type angle)
{
    switch (order)
    {
    case 0x0: // x-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = std::sin(angle * static_cast<m_type>(0.5));
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 0x1: // y-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = 0;
        m_elem[2] = std::sin(angle * static_cast<m_type>(0.5));
        m_elem[3] = 0;
        break;
    case 0x2: // z-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = std::sin(angle * static_cast<m_type>(0.5));
        break;
    default:
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2)
{
    SetElement(order & 0xF, angle1); // Q1
    dtQuaternion<m_type, m_row> Q2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * Q2;
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3)
{
    SetElement(order & 0xF, angle1); // Q1
    dtQuaternion<m_type, m_row> Q2((order >> 4) & 0xF, angle2);
    dtQuaternion<m_type, m_row> Q3((order >> 8) & 0xF, angle3);
    (*this) = (*this) * Q2 * Q3;
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const dtQuaternion &q)
{
    m_elem[0] = q.m_elem[0];
    m_elem[1] = q.m_elem[1];
    m_elem[2] = q.m_elem[2];
    m_elem[3] = q.m_elem[3];
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const uint16_t order, const dtVector3<m_type, 3> &e)
{
    Euler2Quat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const uint16_t order, const dtVector<3, m_type> &e)
{
    Euler2Quat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row>::dtQuaternion(const dtRotation<m_type, 3, 3> &rm)
{
    RotMat2Quat(rm.m_elem);
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetZero()
{
    m_elem[0] = 1;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetFill(const m_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const m_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const m_type w, const m_type x, const m_type y, const m_type z)
{
    m_elem[0] = w;
    m_elem[1] = x;
    m_elem[2] = y;
    m_elem[3] = z;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const uint16_t order, const m_type angle)
{
    switch (order)
    {
    case 0x0: // x-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = std::sin(angle * static_cast<m_type>(0.5));
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 0x1: // y-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = 0;
        m_elem[2] = std::sin(angle * static_cast<m_type>(0.5));
        m_elem[3] = 0;
        break;
    case 0x2: // z-axis
        m_elem[0] = std::cos(angle * static_cast<m_type>(0.5));
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = std::sin(angle * static_cast<m_type>(0.5));
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const uint16_t order, const m_type angle1, const m_type angle2)
{
    SetElement(order & 0xF, angle1); // Q1
    dtQuaternion<m_type, m_row> Q2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * Q2;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3)
{
    SetElement(order & 0xF, angle1); // Q1
    dtQuaternion<m_type, m_row> Q2((order >> 4) & 0xF, angle2);
    dtQuaternion<m_type, m_row> Q3((order >> 8) & 0xF, angle3);
    (*this) = (*this) * Q2 * Q3;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const dtQuaternion &q)
{
    m_elem[0] = q.m_elem[0];
    m_elem[1] = q.m_elem[1];
    m_elem[2] = q.m_elem[2];
    m_elem[3] = q.m_elem[3];
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const uint16_t order, const dtVector3<m_type, 3> &e)
{
    Euler2Quat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const uint16_t order, const dtVector<3, m_type> &e)
{
    Euler2Quat(order, e.m_elem);
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetElement(const dtRotation<m_type, 3, 3> &rm)
{
    RotMat2Quat(rm.m_elem);
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetSwap(const uint16_t i, const uint16_t j)
{
    m_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::SetNormalize()
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    m_elem[0] /= norm;
    m_elem[1] /= norm;
    m_elem[2] /= norm;
    m_elem[3] /= norm;
}

template <typename m_type, uint16_t m_row>
inline const m_type *const dtQuaternion<m_type, m_row>::GetElementsAddr() const
{
    return m_elem;
}

template <typename m_type, uint16_t m_row>
inline m_type dtQuaternion<m_type, m_row>::GetNorm() const
{
    return std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtQuaternion<m_type, m_row>::GetSqNorm() const
{
    return (
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline m_type dtQuaternion<m_type, m_row>::GetSum() const
{
    return (
        m_elem[0] +
        m_elem[1] +
        m_elem[2] +
        m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::GetNormalized() const
{
    m_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<m_type>::epsilon())
        norm = std::numeric_limits<m_type>::epsilon();

    return dtQuaternion(
        m_elem[0] / norm,
        m_elem[1] / norm,
        m_elem[2] / norm,
        m_elem[3] / norm);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::GetConj() const
{
    return dtQuaternion(
        m_elem[0],
        -m_elem[1],
        -m_elem[2],
        -m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtQuaternion<m_type, m_row>::GetEulerAngles(const uint16_t order) const
{
    /* Tait-Bryan angles */
    m_type vec[3];
    m_type pivot;

    // 0:x, 1:y, 2:z
    // order = 0x012 -> zyx, 0x210 -> xyz, 0x102 -> zxy, inverse order!
    uint32_t o1 = (order & 0xF) + 1;
    uint32_t o2 = ((order >> 4) & 0xF) + 1;
    uint32_t o3 = ((order >> 8) & 0xF) + 1;
    int sign = ((o1 + 1) == o2) ? 1 : -1;

    pivot = sign * 2 * (m_elem[o1] * m_elem[o3] + sign * m_elem[0] * m_elem[o2]);
    vec[1] = std::asin(pivot);

    if ((1 - std::fabs(pivot)) <= std::numeric_limits<m_type>::epsilon())
    {
        vec[0] = std::atan2(
            sign * 2 * (m_elem[o3] * m_elem[o2] + sign * m_elem[0] * m_elem[o1]),
            1 - 2 * (m_elem[o1] * m_elem[o1] + m_elem[o3] * m_elem[o3]));
        vec[2] = 0;
    }
    else
    {
        vec[0] = std::atan2(
            -sign * 2 * (m_elem[o3] * m_elem[o2] - sign * m_elem[0] * m_elem[o1]),
            1 - 2 * (m_elem[o1] * m_elem[o1] + m_elem[o2] * m_elem[o2]));
        vec[2] = std::atan2(
            -sign * 2 * (m_elem[o1] * m_elem[o2] - sign * m_elem[0] * m_elem[o3]),
            1 - 2 * (m_elem[o3] * m_elem[o3] + m_elem[o2] * m_elem[o2]));
    }

    return dtVector3<m_type, 3>(vec);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtQuaternion<m_type, m_row>::GetOriErr(const dtQuaternion &q) const
{
    ////Ref: S.-K. Kim's "Concurrent control of position / orientation of a redundant manipulator based on virtual springdamper hypothesis"
    ///* Orientation error (X - Xd) in rotation matrix */
    //// Re = Rd.Transpose() * R : Re is error R, Rd is desired R, R is current R
    ///* Orientation error (X - Xd) in quaternion */
    //// qe = qd.Conj() * q : qe is error q, qd is desired q, q is current q
    //
    ///* Step1. Get Orientation error in quaternion */
    //// qe = qd.Conj() * q
    //// desired q is this quaternion, current q is argument q
    // m_type qe[m_row];
    // qe[0] = m_elem[0] * q.m_elem[0] + m_elem[1] * q.m_elem[1] + m_elem[2] * q.m_elem[2] + m_elem[3] * q.m_elem[3];
    // qe[1] = m_elem[0] * q.m_elem[1] - m_elem[1] * q.m_elem[0] - m_elem[2] * q.m_elem[3] + m_elem[3] * q.m_elem[2];
    // qe[2] = m_elem[0] * q.m_elem[2] + m_elem[1] * q.m_elem[3] - m_elem[2] * q.m_elem[0] - m_elem[3] * q.m_elem[1];
    // qe[3] = m_elem[0] * q.m_elem[3] - m_elem[1] * q.m_elem[2] + m_elem[2] * q.m_elem[1] - m_elem[3] * q.m_elem[0];
    //
    ///* Step2. Orientatin error for geometric jacobian from quaternion */
    //// qe = (1/2)*T.Transpose()*qe(1:3) or (1/2)*qe(0)*qe(0:3)
    //// T = qe(0)*I3 + qe(1:3).GetSkew() : qe(1:3) is vector that consist of qe(1), qe(2) and qe(3). I3 is identity matrix
    //// qe means 'orientation - desired orentation', so qe needed '-' operator to be 'desired orientatin - orientation'
    // return dtVector3<m_type, 3>(-0.5f * qe[0] * qe[1], -0.5f * qe[0] * qe[2], -0.5f * qe[0] * qe[3]);

    /* Old version */
    // dtQuaternion<m_type> quatErr = (*this) / q;
    dtQuaternion<m_type, 4> quatErr; // desired q conjugate * q
    quatErr.SetElement(
        m_elem[0] * q.m_elem[0] + m_elem[1] * q.m_elem[1] + m_elem[2] * q.m_elem[2] + m_elem[3] * q.m_elem[3],
        -m_elem[0] * q.m_elem[1] + m_elem[1] * q.m_elem[0] - m_elem[2] * q.m_elem[3] + m_elem[3] * q.m_elem[2],
        -m_elem[0] * q.m_elem[2] + m_elem[1] * q.m_elem[3] + m_elem[2] * q.m_elem[0] - m_elem[3] * q.m_elem[1],
        -m_elem[0] * q.m_elem[3] - m_elem[1] * q.m_elem[2] + m_elem[2] * q.m_elem[1] + m_elem[3] * q.m_elem[0]);

    dtMatrix3<m_type, 3, 3> quatErr0Mat(
        quatErr.m_elem[0], 0, 0,
        0, quatErr.m_elem[0], 0,
        0, 0, quatErr.m_elem[0]);

    dtVector3<m_type, 3> quatErrEps(quatErr.m_elem[1], quatErr.m_elem[2], quatErr.m_elem[3]);

    dtMatrix3<m_type, 3, 3> U = (quatErr0Mat + quatErrEps.GetSkew()) * 0.5;

    return dtVector3<m_type>(U.Transpose() * quatErrEps);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::exp() const
{
    /* Exponential of general quaternions */
    // e^(q) = e^(qw) * e^(qv) = e^(qw) * e^(u*th)
    // e^(q) = exp(qw) * [cos(|qv|) u * sin(|qv|)]T
    //       = exp(qw) * [cos(|qv|) sin(|qv|)*qx/|qv| sin(|qv|)*qy/|qv| sin(|qv|)*qz/|qv|]T
    // where q = [qw qx qy qz]T, qv = [qx qy qz]T, u is unit vector, u*th = qv
    //
    // q = exp(u*phi/2)
    // u * phi is rotation vector

    m_type e_qw = std::exp(m_elem[0]);
    m_type norm_qv = std::sqrt(m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type alpha = e_qw * sin(norm_qv) / norm_qv;

    return dtQuaternion(
        e_qw * std::cos(norm_qv),
        m_elem[1] * alpha,
        m_elem[2] * alpha,
        m_elem[3] * alpha);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::log() const
{
    /* Logarithm of general quaternions */
    // log(q) = ln(q) = [ln(|q|) th*qx/|qv| th*qy/|qv| th*qz/|qv|]T = [log(|q|) u*th]
    // if |q| = 1, log(q) = [0 th*qx/|qv| th*qy/|qv| th*qz/|qv|]T = [0 u*th]
    // where q = [qw qx qy qz]T, qv = [qx qy qz]T, th = arctan2(|qv|, qw)

    m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type norm_qv = std::sqrt(m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type alpha;

    if (norm_qv > std::numeric_limits<m_type>::epsilon())
        alpha = std::atan2(norm_qv, m_elem[0]) / norm_qv; // th / norm_qv
    else
        alpha = 0; // singular

    return dtQuaternion(
        std::log(norm_q),
        m_elem[1] * alpha,
        m_elem[2] * alpha,
        m_elem[3] * alpha);
}

template <typename m_type, uint16_t m_row>
inline dtVector3<m_type, 3> dtQuaternion<m_type, m_row>::Log() const
{
    /* Capitalized Logarithm of quaternions */
    // S3 -> R3
    // Log(q) = u*phi = u*2*th, when q is unit quternion

    // m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    // m_type unit_q[4] = { m_elem[0] / norm_q, m_elem[1] / norm_q, m_elem[2] / norm_q, m_elem[3] / norm_q };
    // m_type norm_qv = std::sqrt(unit_q[1] * unit_q[1] + unit_q[2] * unit_q[2] + unit_q[3] * unit_q[3]);
    // m_type alpha = 2 * std::atan2(norm_qv, unit_q[0]) / norm_qv; // phi / norm_qv
    // return dtVector3<m_type, 3>(unit_q[1] * alpha, unit_q[2] * alpha, unit_q[3] * alpha);

    m_type norm_qv = std::sqrt(m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type alpha;

    if (norm_qv > std::numeric_limits<m_type>::epsilon())
        alpha = 2 * std::atan2(norm_qv, m_elem[0]) / norm_qv; // phi / norm_qv
    else
        alpha = 0; // singular

    return dtVector3<m_type, 3>(m_elem[1] * alpha, m_elem[2] * alpha, m_elem[3] * alpha);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::ode(m_type wx, m_type wy, m_type wz) const
{
    /* Ordinary Differential Equation (ODE) */
    // dq/dt = (q * Wl) / 2 = (Wg * q) / 2
    // where Wl is local angular velocity, Wg is global angular velocity
    m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type unit_q[4] = {m_elem[0] / norm_q, m_elem[1] / norm_q, m_elem[2] / norm_q, m_elem[3] / norm_q};

    return dtQuaternion(
        (-unit_q[1] * wx - unit_q[2] * wy - unit_q[3] * wz) * 0.5,
        (unit_q[0] * wx + unit_q[2] * wz - unit_q[3] * wy) * 0.5,
        (unit_q[0] * wy - unit_q[1] * wz + unit_q[3] * wx) * 0.5,
        (unit_q[0] * wz + unit_q[1] * wy - unit_q[2] * wx) * 0.5);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::ode(m_type *w) const
{
    /* Ordinary Differential Equation (ODE) */
    // dq/dt = (q * Wl) / 2 = (Wg * q) / 2
    // where Wl is local angular velocity, Wg is global angular velocity
    m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type unit_q[4] = {m_elem[0] / norm_q, m_elem[1] / norm_q, m_elem[2] / norm_q, m_elem[3] / norm_q};

    return dtQuaternion(
        (-unit_q[1] * w[0] - unit_q[2] * w[1] - unit_q[3] * w[2]) * 0.5,
        (unit_q[0] * w[0] + unit_q[2] * w[2] - unit_q[3] * w[1]) * 0.5,
        (unit_q[0] * w[1] - unit_q[1] * w[2] + unit_q[3] * w[0]) * 0.5,
        (unit_q[0] * w[2] + unit_q[1] * w[1] - unit_q[2] * w[0]) * 0.5);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::ode(dtVector3<m_type, 3> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // dq/dt = (q * Wl) / 2 = (Wg * q) / 2
    // where Wl is local angular velocity, Wg is global angular velocity
    m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type unit_q[4] = {m_elem[0] / norm_q, m_elem[1] / norm_q, m_elem[2] / norm_q, m_elem[3] / norm_q};

    return dtQuaternion(
        (-unit_q[1] * w.m_elem[0] - unit_q[2] * w.m_elem[1] - unit_q[3] * w.m_elem[2]) * 0.5,
        (unit_q[0] * w.m_elem[0] + unit_q[2] * w.m_elem[2] - unit_q[3] * w.m_elem[1]) * 0.5,
        (unit_q[0] * w.m_elem[1] - unit_q[1] * w.m_elem[2] + unit_q[3] * w.m_elem[0]) * 0.5,
        (unit_q[0] * w.m_elem[2] + unit_q[1] * w.m_elem[1] - unit_q[2] * w.m_elem[0]) * 0.5);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::ode(dtVector<3, m_type> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // dq/dt = (q * Wl) / 2 = (Wg * q) / 2
    // where Wl is local angular velocity, Wg is global angular velocity
    m_type norm_q = std::sqrt(m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3]);
    m_type unit_q[4] = {m_elem[0] / norm_q, m_elem[1] / norm_q, m_elem[2] / norm_q, m_elem[3] / norm_q};

    return dtQuaternion(
        (-unit_q[1] * w.m_elem[0] - unit_q[2] * w.m_elem[1] - unit_q[3] * w.m_elem[2]) * 0.5,
        (unit_q[0] * w.m_elem[0] + unit_q[2] * w.m_elem[2] - unit_q[3] * w.m_elem[1]) * 0.5,
        (unit_q[0] * w.m_elem[1] - unit_q[1] * w.m_elem[2] + unit_q[3] * w.m_elem[0]) * 0.5,
        (unit_q[0] * w.m_elem[2] + unit_q[1] * w.m_elem[1] - unit_q[2] * w.m_elem[0]) * 0.5);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::Inv() const
{
    /* Inverse of quaternion */
    // q^(-1) = q.Conj() / (q.GetNorm())^2
    // if |q| = 1, q^(-1) = q.Conj() : unit quaternion

    m_type norm2 = m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] + m_elem[3] * m_elem[3];

    return dtQuaternion(
        m_elem[0] / norm2,
        -m_elem[1] / norm2,
        -m_elem[2] / norm2,
        -m_elem[3] / norm2);
}

/* Assignment operators */
template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> &dtQuaternion<m_type, m_row>::operator=(const dtQuaternion &q)
{
    m_elem[0] = q.m_elem[0];
    m_elem[1] = q.m_elem[1];
    m_elem[2] = q.m_elem[2];
    m_elem[3] = q.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> &dtQuaternion<m_type, m_row>::operator+=(const dtQuaternion &q)
{
    m_elem[0] += q.m_elem[0];
    m_elem[1] += q.m_elem[1];
    m_elem[2] += q.m_elem[2];
    m_elem[3] += q.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> &dtQuaternion<m_type, m_row>::operator-=(const dtQuaternion &q)
{
    m_elem[0] -= q.m_elem[0];
    m_elem[1] -= q.m_elem[1];
    m_elem[2] -= q.m_elem[2];
    m_elem[3] -= q.m_elem[3];

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> &dtQuaternion<m_type, m_row>::operator*=(const m_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> &dtQuaternion<m_type, m_row>::operator/=(const m_type s)
{
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    m_elem[0] /= scalar;
    m_elem[1] /= scalar;
    m_elem[2] /= scalar;
    m_elem[3] /= scalar;

    return (*this);
}

template <typename m_type, uint16_t m_row>
inline dtCommaInit<m_row, m_type> dtQuaternion<m_type, m_row>::operator<<(const m_type s)
{
    m_elem[0] = s;
    return dtCommaInit<m_row, m_type>(m_elem);
}

/* Arithmetic operators */
template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator-() const
{
    return dtQuaternion(
        -m_elem[0],
        -m_elem[1],
        -m_elem[2],
        -m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator+(const dtQuaternion &q) const
{
    return dtQuaternion(
        m_elem[0] + q.m_elem[0],
        m_elem[1] + q.m_elem[1],
        m_elem[2] + q.m_elem[2],
        m_elem[3] + q.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator-(const dtQuaternion &q) const
{
    return dtQuaternion(
        m_elem[0] - q.m_elem[0],
        m_elem[1] - q.m_elem[1],
        m_elem[2] - q.m_elem[2],
        m_elem[3] - q.m_elem[3]);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator*(const m_type s) const
{
    return dtQuaternion(
        m_elem[0] * s,
        m_elem[1] * s,
        m_elem[2] * s,
        m_elem[3] * s);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator/(const m_type s) const
{
    m_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<m_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<m_type>::epsilon();
        else
            scalar = std::numeric_limits<m_type>::epsilon();
    }

    return dtQuaternion(
        m_elem[0] / scalar,
        m_elem[1] / scalar,
        m_elem[2] / scalar,
        m_elem[3] / scalar);
}

template <typename m_type, uint16_t m_row>
inline dtQuaternion<m_type, m_row> dtQuaternion<m_type, m_row>::operator*(const dtQuaternion &q) const
{
    return dtQuaternion<m_type, m_row>(
        m_elem[0] * q.m_elem[0] - m_elem[1] * q.m_elem[1] - m_elem[2] * q.m_elem[2] - m_elem[3] * q.m_elem[3],
        m_elem[0] * q.m_elem[1] + m_elem[1] * q.m_elem[0] + m_elem[2] * q.m_elem[3] - m_elem[3] * q.m_elem[2],
        m_elem[0] * q.m_elem[2] - m_elem[1] * q.m_elem[3] + m_elem[2] * q.m_elem[0] + m_elem[3] * q.m_elem[1],
        m_elem[0] * q.m_elem[3] + m_elem[1] * q.m_elem[2] - m_elem[2] * q.m_elem[1] + m_elem[3] * q.m_elem[0]);
}

/* Comparison operators */
template <typename m_type, uint16_t m_row>
inline bool dtQuaternion<m_type, m_row>::operator==(const dtQuaternion &q) const
{
    // q and -q are same, because quaternion double covered
    if (SGN(m_elem[0]) == SGN(q.m_elem[0]))
    {
        if (m_elem[0] - q.m_elem[0] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[1] - q.m_elem[1] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[2] - q.m_elem[2] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[3] - q.m_elem[3] > std::numeric_limits<m_type>::epsilon())
            return false;
    }
    else
    {
        if (m_elem[0] + q.m_elem[0] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[1] + q.m_elem[1] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[2] + q.m_elem[2] > std::numeric_limits<m_type>::epsilon())
            return false;
        else if (m_elem[3] + q.m_elem[3] > std::numeric_limits<m_type>::epsilon())
            return false;
    }

    return true;
}

template <typename m_type, uint16_t m_row>
inline bool dtQuaternion<m_type, m_row>::operator!=(const dtQuaternion &q) const
{
    // q and -q are same, because quaternion double covered
    if (SGN(m_elem[0]) == SGN(q.m_elem[0]))
    {
        if (m_elem[0] - q.m_elem[0] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[1] - q.m_elem[1] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[2] - q.m_elem[2] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[3] - q.m_elem[3] > std::numeric_limits<m_type>::epsilon())
            return true;
    }
    else
    {
        if (m_elem[0] + q.m_elem[0] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[1] + q.m_elem[1] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[2] + q.m_elem[2] > std::numeric_limits<m_type>::epsilon())
            return true;
        else if (m_elem[3] + q.m_elem[3] > std::numeric_limits<m_type>::epsilon())
            return true;
    }

    return false;
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        Serial.printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        printf("%7.3f\n", (m_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Private Member Function ------------------------------------------------//
template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::Euler2Quat(const uint16_t order, const m_type *e)
{
    m_type s_ps = std::sin(e[0] * static_cast<m_type>(0.5)); // sin(psi)
    m_type c_ps = std::cos(e[0] * static_cast<m_type>(0.5)); // cos(psi)
    m_type s_th = std::sin(e[1] * static_cast<m_type>(0.5)); // sin(the)
    m_type c_th = std::cos(e[1] * static_cast<m_type>(0.5)); // cos(the)
    m_type s_ph = std::sin(e[2] * static_cast<m_type>(0.5)); // sin(phi)
    m_type c_ph = std::cos(e[2] * static_cast<m_type>(0.5)); // cos(phi)

    /* Only Tait?Bryan angles */
    switch (order)
    {
    case 0x120: // xzy
        m_elem[0] = s_ps * s_th * s_ph + c_ps * c_th * c_ph;
        m_elem[1] = s_ps * c_th * c_ph - c_ps * s_th * s_ph;
        m_elem[2] = c_ps * c_th * s_ph - s_ps * s_th * c_ph;
        m_elem[3] = s_ps * c_th * s_ph + c_ps * s_th * c_ph;
        break;
    case 0x210: // xyz
        m_elem[0] = c_ps * c_th * c_ph - s_ps * s_th * s_ph;
        m_elem[1] = c_ps * s_th * s_ph + s_ps * c_th * c_ph;
        m_elem[2] = c_ps * s_th * c_ph - s_ps * c_th * s_ph;
        m_elem[3] = c_ps * c_th * s_ph + s_ps * s_th * c_ph;
        break;
    case 0x201: // yxz
        m_elem[0] = s_ps * s_th * s_ph + c_ps * c_th * c_ph;
        m_elem[1] = s_ps * c_th * s_ph + c_ps * s_th * c_ph;
        m_elem[2] = s_ps * c_th * c_ph - c_ps * s_th * s_ph;
        m_elem[3] = c_ps * c_th * s_ph - s_ps * s_th * c_ph;
        break;
    case 0x021: // yzx
        m_elem[0] = c_ps * c_th * c_ph - s_ps * s_th * s_ph;
        m_elem[1] = c_ps * c_th * s_ph + s_ps * s_th * c_ph;
        m_elem[2] = c_ps * s_th * s_ph + s_ps * c_th * c_ph;
        m_elem[3] = c_ps * s_th * c_ph - s_ps * c_th * s_ph;
        break;
    case 0x012: // zyx
        m_elem[0] = s_ps * s_th * s_ph + c_ps * c_th * c_ph;
        m_elem[1] = c_ps * c_th * s_ph - s_ps * s_th * c_ph;
        m_elem[2] = s_ps * c_th * s_ph + c_ps * s_th * c_ph;
        m_elem[3] = s_ps * c_th * c_ph - c_ps * s_th * s_ph;
        break;
    case 0x102: // zxy
        m_elem[0] = c_ps * c_th * c_ph - s_ps * s_th * s_ph;
        m_elem[1] = c_ps * s_th * c_ph - s_ps * c_th * s_ph;
        m_elem[2] = c_ps * c_th * s_ph + s_ps * s_th * c_ph;
        m_elem[3] = c_ps * s_th * s_ph + s_ps * c_th * c_ph;
        break;
    }
}

template <typename m_type, uint16_t m_row>
inline void dtQuaternion<m_type, m_row>::RotMat2Quat(const m_type *rm)
{
    // Get squared elements
    m_elem[0] = (1 + rm[0] + rm[4] + rm[8]) * static_cast<m_type>(0.25); // w^2
    m_elem[1] = (1 + rm[0] - rm[4] - rm[8]) * static_cast<m_type>(0.25); // x^2
    m_elem[2] = (1 - rm[0] + rm[4] - rm[8]) * static_cast<m_type>(0.25); // y^2
    m_elem[3] = (1 - rm[0] - rm[4] + rm[8]) * static_cast<m_type>(0.25); // z^2

    // Get element value but this value is always positive
    m_elem[0] = std::sqrt(m_elem[0]); // |w|
    m_elem[1] = std::sqrt(m_elem[1]); // |x|
    m_elem[2] = std::sqrt(m_elem[2]); // |y|
    m_elem[3] = std::sqrt(m_elem[3]); // |z|

    // Choose the sign of the quaternion's element
    if (m_elem[0] <= std::numeric_limits<m_type>::epsilon())
    {
        if (m_elem[1] <= std::numeric_limits<m_type>::epsilon())
        {
            if (m_elem[2] <= std::numeric_limits<m_type>::epsilon())
            { /* w == 0 && x == 0 && y == 0*/
                m_elem[0] = 0;
                m_elem[1] = 0;
                m_elem[2] = 0;
            }
            else
            { /* w == 0 && x == 0 */
                m_elem[0] = 0;
                m_elem[1] = 0;
                m_elem[3] *= SGN(rm[7] + rm[5]); // y*z = (m21 + m12) / 4
            }
        }
        else
        { /* w == 0 */
            m_elem[0] = 0;
            m_elem[2] *= SGN(rm[3] + rm[1]); // x*y = (m10 + m01) / 4
            m_elem[3] *= SGN(rm[2] + rm[6]); // x*z = (m02 + m20) / 4
        }
    }
    else
    {
        m_elem[1] *= SGN(rm[7] - rm[5]); // w*x = (m21 - m12) / 4
        m_elem[2] *= SGN(rm[2] - rm[6]); // w*y = (m02 - m20) / 4
        m_elem[3] *= SGN(rm[3] - rm[1]); // w*z = (m10 - m01) / 4
    }
}

//-- Template Function ------------------------------------------------------//
// scalar * quaternion
template <typename type, uint16_t row>
inline dtQuaternion<type, row> operator*(const type s, const dtQuaternion<type, row> &q)
{
    return dtQuaternion<type, row>(
        q.m_elem[0] * s,
        q.m_elem[1] * s,
        q.m_elem[2] * s,
        q.m_elem[3] * s);
}

typedef dtQuaternion<> dtQuat;

} // namespace dtMath

#endif // DTMATH_DTQUATERNION_TPP_