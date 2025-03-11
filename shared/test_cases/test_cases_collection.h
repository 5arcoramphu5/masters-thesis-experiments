#pragma once

#include "test_case.h"

class TestCasesCollection
{
    private:
        int maxDerivative;

        void generate_diagonal_matrix_test();
        void generate_PCR3BP_test();

    public:
        TestCase diagonal_matrix;
        TestCase PCR3BP;

        TestCasesCollection(int maxDerivative) :  maxDerivative(maxDerivative)
        {
            generate_diagonal_matrix_test();
            generate_PCR3BP_test();
        }
};