/*!
\file       dtContiRiccatiEq.h
\brief      dtMath, The continuous-time algebraic Riccati Equation Class, cmake version
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\see        https://github.com/RobotLocomotion/drake/tree/master/math
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCARE_TPP_
#define DTMATH_DTCARE_TPP_

#include "dtCARE.h"

namespace dtMath
{

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline dtCARE<m_dimN, m_dimM, m_type>::dtCARE()
{
    m_mEye.SetIdentity();
    m_maxIteration = 100;
    m_tolerance = std::numeric_limits<m_type>::epsilon();
    // m_tolerance = 1e-9;
}

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline void dtCARE<m_dimN, m_dimM, m_type>::SetSolveCriteria(m_type tolerance, uint16_t maxIteration)
{
    m_tolerance = tolerance;
    m_maxIteration = maxIteration;
}

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline dtMatrix<m_dimN, m_dimN, m_type> dtCARE<m_dimN, m_dimM, m_type>::Solve(
    const dtMatrix<m_dimN, m_dimN, m_type> &A,
    const dtMatrix<m_dimN, m_dimM, m_type> &B,
    const dtMatrix<m_dimN, m_dimN, m_type> &Q,
    const dtMatrix<m_dimM, m_dimM, m_type> &R)
{
    // Is Q Symmetric matrix?
    if (Q != Q.Transpose())
        return dtMatrix<m_dimN, m_dimN, m_type>();

    m_type relative_norm;
    uint16_t iteration = 0;

    // m_mH.SetBlock(0, 0, A);
    // m_mH.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    // m_mH.SetBlock(m_dimN, 0, Q);
    // m_mH.SetBlock(m_dimN, m_dimN, -A.Transpose());
    //
    // m_mZ = m_mH;

    m_mZ.SetBlock(0, 0, A);
    m_mZ.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    m_mZ.SetBlock(m_dimN, 0, Q);
    m_mZ.SetBlock(m_dimN, m_dimN, -A.Transpose());

    do
    {
        m_mZold = m_mZ;
        // R. Byers. Solving the algebraic Riccati equation with the matrix sign
        // function. Linear Algebra Appl., 85:267?279, 1987
        // Added determinant scaling to improve convergence (converges in rough half
        // the iterations with this)
        m_type ck = std::pow(std::abs(m_mZ.Determinant()), m_power);
        m_mZ *= ck;
        m_mZ = m_mZ - (m_mZ - m_mZ.Inv()) * 0.5;
        relative_norm = (m_mZ - m_mZold).GetNorm();
        iteration++;
    } while ((iteration < m_maxIteration) && (relative_norm > m_tolerance));

    m_mZ.GetBlock(0, 0, m_mW11);
    m_mZ.GetBlock(0, m_dimN, m_mW12);
    m_mZ.GetBlock(m_dimN, 0, m_mW21);
    m_mZ.GetBlock(m_dimN, m_dimN, m_mW22);

    m_mLhs.SetBlock(0, 0, m_mW12);
    m_mLhs.SetBlock(m_dimN, 0, m_mW22 + m_mEye);

    m_mRhs.SetBlock(0, 0, m_mW11 + m_mEye);
    m_mRhs.SetBlock(m_dimN, 0, m_mW21);

    return m_mLhs.SVD().Solve(m_mRhs);
}

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline int8_t dtCARE<m_dimN, m_dimM, m_type>::Solve(
    const dtMatrix<m_dimN, m_dimN, m_type> &A,
    const dtMatrix<m_dimN, m_dimM, m_type> &B,
    const dtMatrix<m_dimN, m_dimN, m_type> &Q,
    const dtMatrix<m_dimM, m_dimM, m_type> &R,
    dtMatrix<m_dimN, m_dimN, m_type> &X)
{
    // Is Q Symmetric matrix?
    if (Q != Q.Transpose())
        return -1;

    m_type relative_norm;
    uint16_t iteration = 0;

    // m_mH.SetBlock(0, 0, A);
    // m_mH.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    // m_mH.SetBlock(m_dimN, 0, Q);
    // m_mH.SetBlock(m_dimN, m_dimN, -A.Transpose());
    //
    // m_mZ = m_mH;

    m_mZ.SetBlock(0, 0, A);
    m_mZ.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    m_mZ.SetBlock(m_dimN, 0, Q);
    m_mZ.SetBlock(m_dimN, m_dimN, -A.Transpose());

    do
    {
        m_mZold = m_mZ;
        // R. Byers. Solving the algebraic Riccati equation with the matrix sign
        // function. Linear Algebra Appl., 85:267?279, 1987
        // Added determinant scaling to improve convergence (converges in rough half
        // the iterations with this)
        m_type ck = std::pow(std::abs(m_mZ.Determinant()), m_power);
        m_mZ *= ck;
        m_mZ = m_mZ - (m_mZ - m_mZ.Inv()) * 0.5;
        relative_norm = (m_mZ - m_mZold).GetNorm();
        iteration++;
    } while ((iteration < m_maxIteration) && (relative_norm > m_tolerance));

    m_mZ.GetBlock(0, 0, m_mW11);
    m_mZ.GetBlock(0, m_dimN, m_mW12);
    m_mZ.GetBlock(m_dimN, 0, m_mW21);
    m_mZ.GetBlock(m_dimN, m_dimN, m_mW22);

    m_mLhs.SetBlock(0, 0, m_mW12);
    m_mLhs.SetBlock(m_dimN, 0, m_mW22 + m_mEye);

    m_mRhs.SetBlock(0, 0, m_mW11 + m_mEye);
    m_mRhs.SetBlock(m_dimN, 0, m_mW21);

    return m_mLhs.SVD().Solve(m_mRhs, X);
}

} // namespace dtMath

#endif // DTMATH_DTCARE_TPP_