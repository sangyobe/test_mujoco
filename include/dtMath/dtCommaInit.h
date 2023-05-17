/*!
\file       dtCommaInit.h
\brief      dtMath, Matrix & Vector comma initializer class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history!
*/

#ifndef DTMATH_DTCOMMA_INIT_H_
#define DTMATH_DTCOMMA_INIT_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

namespace dtMath
{

template <uint16_t m_size, typename m_type = float>
class dtCommaInit
{
private:
    m_type *m_elem;
    int m_idx = 1;

public:
    dtCommaInit(m_type *elem) : m_elem(elem) {}

    dtCommaInit &operator,(const m_type s)
    {
        if (m_idx < m_size)
            m_elem[m_idx++] = s;
        return *this;
    }
};

} // namespace dtMath

#include "dtCommaInit.tpp"

#endif // DTMATH_DTCOMMA_INIT_H_