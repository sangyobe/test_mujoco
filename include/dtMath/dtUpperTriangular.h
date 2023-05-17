/*!
\file       dtUpperTriangular.h
\brief      dtMath, Upper triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTUPPER_TRIANGULAR_H_
#define DTMATH_DTUPPER_TRIANGULAR_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dtMath
{

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class dtUpperTriangular
{
public:
    dtUpperTriangular() {}

    // for General Upper Triangular Matrix
    static int8_t Solve(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);
    static dtVector<m_col, m_type> Solve(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t Inverse(dtMatrix<m_row, m_col, m_type> U, dtMatrix<m_row, m_col, m_type> &invU);
    static dtMatrix<m_row, m_col, m_type> Inverse(dtMatrix<m_row, m_col, m_type> U, int8_t *isOk = nullptr);

    // for Unit Upper Triangular Matrix
    static int8_t SolveUnit(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);
    static dtVector<m_col, m_type> SolveUnit(dtMatrix<m_row, m_col, m_type> &U, dtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t InverseUnit(dtMatrix<m_row, m_col, m_type> U, dtMatrix<m_row, m_col, m_type> &invU);
    static dtMatrix<m_row, m_col, m_type> InverseUnit(dtMatrix<m_row, m_col, m_type> U, int8_t *isOk = nullptr);
};

} // namespace dtMath

#include "dtUpperTriangular.tpp"

#endif // DTMATH_DTUPPER_TRIANGULAR_H_
