#include <gtest/gtest.h>
#include <dtMath/dtMath.h>

template <uint16_t m_row, uint16_t m_col, typename m_type>
std::ostream &operator<<(std::ostream &os, dtMath::dtMatrix<m_row, m_col, m_type> &mat)
{
    for (uint16_t i = 0; i < m_row; i++)
    {
        for (uint16_t j = 0; j < m_col; j++)
        {
            std::cout << mat(i, j) << " ";
        }
        std::cout << std::endl;
    }
    return os;
}

template <>
std::ostream &operator<< <3, 3, double>(std::ostream &os, dtMath::dtMatrix<3, 3, double> &mat)
{
    std::cout << "3x3 mat = " << std::endl;
    std::cout << mat(0, 0) << " " << mat(0, 1) << " " << mat(0, 2) << std::endl;
    std::cout << mat(1, 0) << " " << mat(1, 1) << " " << mat(1, 2) << std::endl;
    std::cout << mat(2, 0) << " " << mat(2, 1) << " " << mat(2, 2) << std::endl;
    return os;
}

namespace TestMath
{
    namespace
    {
        TEST(TestMath, Initialize)
        {
            dtMath::dtMatrix<3, 3, double> m1;
            std::cout << m1 << std::endl;
        }
    }
}