#pragma once

#include "../../../normal-forms/source/typedefs.h"
#include "../../../normal-forms/source/Diagonalization/Diagonalization.hpp"

class TestCasesCollection
{
    private:
        int maxDerivative;

        void generate_diagonal_matrix_test();
        void generate_PCR3BP_test();

    public:
        Diagonalization<capd::Complex> diagonal_matrix;
        Diagonalization<capd::Complex> PCR3BP;

        capd::DMap PCR3BP_dmap;

        TestCasesCollection(int maxDerivative); 
};