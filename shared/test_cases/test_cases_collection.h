#pragma once

#include "../../../normal-forms/source/typedefs.h"
#include "../../../normal-forms/source/Diagonalization/Diagonalization.hpp"

template<ArithmeticType Coeff>
struct TestCase
{
    std::string testName;
    Diagonalization<Coeff> diagonalization;

    TestCase(std::string testName, MapFunction f, int noParams, const Vector<Coeff> &p, const Matrix<Coeff> &J, const Matrix<Coeff> &invJ, const Matrix<Coeff> &lambda, int maxDerivative)
     : testName(testName), diagonalization(f, noParams, p, J, invJ, lambda, maxDerivative) {}
     
    TestCase() : testName(""), diagonalization() {}
};

class TestCasesCollection
{
    private:
        int maxDerivative;

        void generate_diagonal_matrix_test();
        void generate_PCR3BP_L4_test();
        void generate_PCR3BP_L1_test();

    public:
        TestCase<capd::Complex> diagonal_matrix;
        TestCase<capd::Complex> PCR3BP_L4;
        TestCase<capd::Complex> PCR3BP_L1;

        capd::LDMap PCR3BP_dmap;

        TestCasesCollection(int maxDerivative); 
};