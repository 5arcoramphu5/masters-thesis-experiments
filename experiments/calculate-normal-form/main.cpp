#include <fstream>
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;

#define ALG_LOGGER Logger<ProgressIndication>

int main(int argc, char* argv[])
{
    if (argc != 2) throw runtime_error("program expects one argument: method degree"); 

    int methodDegree = atoi( argv[1] );
    cout << "calculating normal form of degree " << methodDegree << endl;

    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    string path = "shared/precalculated-normal-forms/";

    NormalFormFinder<ALG_LOGGER> finder(methodDegree);
    auto normalForm = finder.calculatePseudoNormalForm(testCase.diagonalization);

    cout << "serializing normal form..." << endl;
    ofstream file(path + testCase.testName + "_deg" + to_string(methodDegree) + ".txt");
    normalForm.serialize(file);
    file.close();

    return 0;
}