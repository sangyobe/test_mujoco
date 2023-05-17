/*!
\file       dtLowerTriangular.h
\brief      dtMath, Lower triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLOWER_TRIANGULAR_H_
#define DTMATH_DTLOWER_TRIANGULAR_H_

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
class dtLowerTriangular
{
public:
    dtLowerTriangular() {}

    // for General Lower Triangular Matrix
    static int8_t Solve(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);
    static dtVector<m_col, m_type> Solve(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t Inverse(dtMatrix<m_row, m_col, m_type> L, dtMatrix<m_row, m_col, m_type> &invL);
    static dtMatrix<m_row, m_col, m_type> Inverse(dtMatrix<m_row, m_col, m_type> L, int8_t *isOk = nullptr);

    // for Unit Lower Triangular Matrix
    static int8_t SolveUnit(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, dtVector<m_col, m_type> &x);
    static dtVector<m_col, m_type> SolveUnit(dtMatrix<m_row, m_col, m_type> &L, dtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t InverseUnit(dtMatrix<m_row, m_col, m_type> L, dtMatrix<m_row, m_col, m_type> &invL);
    static dtMatrix<m_row, m_col, m_type> InverseUnit(dtMatrix<m_row, m_col, m_type> L, int8_t *isOk = nullptr);
};

} // namespace dtMath

#include "dtLowerTriangular.tpp"

#endif // DTMATH_DTLOWER_TRIANGULAR_H_
