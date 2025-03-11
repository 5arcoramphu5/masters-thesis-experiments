#include <iostream>

#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;

#define MAX_DERIVATIVE 10
#define METHOD_DEGREE 4

#define LOGGER Logger<Diagnostic, SymbolicPolynomialPrinting, 13>

void runTest(const TestCase &testCase)
{
    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, testCase.f, testCase.p, testCase.lambda, testCase.J, testCase.invJ);
    auto result = finder.calculatePseudoNormalForm();
}

int main()
{
    TestCasesCollection testCases(MAX_DERIVATIVE);

    runTest(testCases.diagonal_matrix);
    runTest(testCases.PCR3BP);

    return 0;
}