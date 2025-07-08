#include <iostream>
#include <fstream>

#include "../../shared/calculations/calculations.h"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;

void performTest(LDVector &realOriginalPoint, const PseudoNormalForm &normalForm, TestCasesCollection &testCases)
{
    int solverDeg = 30;
    LDOdeSolver solver(testCases.PCR3BP_dmap, solverDeg);
    LDTimeMap timeMap(solver);

    double finalTime = 0.5;

    LDVector L4({0, 0.866025403784438646763723170753, 0, 0});

    double dt = 0.1;
    for(double t = 0; t <= finalTime; t += dt)
    {
        cout << "-----------------------------------------------" << endl;
        cout << "t: " << t << endl;

        LDMatrix intDer;
        LDVector point(realOriginalPoint);
        point = timeMap(t, point, intDer);

        LDVector point2(realOriginalPoint);
        LDMatrix nfDer = normalFormOriginalSolutionDerivative(t, point2, normalForm, testCases.PCR3BP_L4.diagonalization);

        cout <<  "diff:\n" << (intDer - nfDer) << endl;

        double maxAbs = 0;
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
                double absDiff = abs(intDer[i][j] - nfDer[i][j]);
                if(absDiff > maxAbs)
                    maxAbs = absDiff;
            }
        cout << maxAbs << endl;
    }
}

int main()
{
    int methodDegree = 5;
    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    cout << setprecision(20) << fixed;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();

    LDVector realOriginalPoint({0, 0.9, 0, 0});
    performTest(realOriginalPoint, normalForm, testCases);

    return 0;
}