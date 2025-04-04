#include <fstream>
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"
#include "../../shared/test_cases/test_cases_collection.h"
#include "../../shared/typedefs.h"

using namespace std;

#define ALG_LOGGER Logger<ProgressIndication>

int main()
{
    int methodDegree = 15;
    TestCasesCollection testCases(16);
    auto &testCase = testCases.PCR3BP_L4;

    string path = "experiments/comparing-trajectories-test/results/";

    NormalFormFinder<ALG_LOGGER> finder(methodDegree,  testCase.diagonalization);
    auto normalForm = finder.calculatePseudoNormalForm();

    ofstream fileN(path + testCase.testName + "_N_deg" + to_string(methodDegree) + ".txt");
    normalForm.getN().serialize(fileN);
    fileN.close();

    ofstream filePhi(path + testCase.testName + "_Phi_deg" + to_string(methodDegree) + ".txt");
    normalForm.getPhi().serialize(filePhi);
    filePhi.close();

    return 0;
}