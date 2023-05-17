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

#ifndef DTMATH_DTCARE_H_
#define DTMATH_DTCARE_H_

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type> class dtMatrix;

// Computes the unique stabilizing solution X to the continuous-time algebraic
// Riccati equation:
//
// equation: X A + A'X - X B R^{-1} B' X + Q = 0
// where: A(n x n), B(n x m), Q(n x n), R(m x m)
//
// @throws std::runtime_error if R is not positive definite.
//
// Based on the Matrix Sign Function method outlined in this paper:
// http://www.engr.iupui.edu/~skoskie/ECE684/Riccati_algorithms.pdf
//

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type = float>
class dtCARE
{
private:
    // dtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mH;
    dtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mZ;
    dtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mZold;

    dtMatrix<m_dimN, m_dimN, m_type> m_mW11;
    dtMatrix<m_dimN, m_dimN, m_type> m_mW12;
    dtMatrix<m_dimN, m_dimN, m_type> m_mW21;
    dtMatrix<m_dimN, m_dimN, m_type> m_mW22;

    dtMatrix<2 * m_dimN, m_dimN> m_mLhs;
    dtMatrix<2 * m_dimN, m_dimN> m_mRhs;
    dtMatrix<m_dimN, m_dimN> m_mEye;

    uint16_t m_maxIteration;
    m_type m_tolerance;
    const m_type m_power = (m_type)(-1) / (2 * m_dimN);

public:
    dtCARE();
    ~dtCARE() {}
    void SetSolveCriteria(m_type tolerance, uint16_t maxIteration);
    dtMatrix<m_dimN, m_dimN, m_type> Solve(
        const dtMatrix<m_dimN, m_dimN, m_type> &A,
        const dtMatrix<m_dimN, m_dimM, m_type> &B,
        const dtMatrix<m_dimN, m_dimN, m_type> &Q,
        const dtMatrix<m_dimM, m_dimM, m_type> &R);

    int8_t Solve(
        const dtMatrix<m_dimN, m_dimN, m_type> &A,
        const dtMatrix<m_dimN, m_dimM, m_type> &B,
        const dtMatrix<m_dimN, m_dimN, m_type> &Q,
        const dtMatrix<m_dimM, m_dimM, m_type> &R,
        dtMatrix<m_dimN, m_dimN, m_type> &X);
};

} // namespace dtMath

#include "dtCARE.tpp"

#endif // DTMATH_DTCARE_H_