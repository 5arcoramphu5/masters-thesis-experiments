#include <iostream>
#include <fstream>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"
#include "../../shared/typedefs.h"

using namespace std;
using namespace capd;

#define ALG_LOGGER Logger<ProgressIndication>

#define INIT_TIME 0
#define FINAL_TIME 1

DTimeMap::SolutionCurve integrateSolution(DMap &map, DVector &point, int order)
{
    DOdeSolver solver(map, order);
    DTimeMap timeMap(solver);

    DTimeMap::SolutionCurve solution(INIT_TIME);
    timeMap(FINAL_TIME, point, solution);
    return solution;
}

int main()
{
    int methodDegree = 15;
    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    NormalFormFinder<ALG_LOGGER> finder(methodDegree);
    auto normalForm = finder.calculatePseudoNormalForm(testCase.diagonalization);

    Complex lambda1(0.0000001, 0), lambda2(-0.0000001, 0); // the closer to 0 the closer to p transformedPoint will be
    CVector point({ lambda1, lambda2, conj(lambda1), conj(lambda2) });
    auto transformedPoint =  testCase.diagonalization.toOriginal(normalForm.getPhi()(point));
    cout << "transformed point: " << transformedPoint << "\n\n";

    DVector realTransformedPoint(4);
    for(int i = 0; i < 4; ++i)
        realTransformedPoint[i] = transformedPoint[i].real();

    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, realTransformedPoint, 30);

    int N = 10;
    double dt = double(FINAL_TIME - INIT_TIME) / N;
    for(double t = INIT_TIME; t <= FINAL_TIME; t += dt)
    {
        auto normalFormSol = normalForm.solution(t, point);
        normalFormSol = testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormSol));

        cout << "t: " << t << " ";
        cout << "diff: ";
        for(int j = 0; j < 4; ++j)
            cout << (intSolution(t)[j] - normalFormSol[j]) << " ";
        cout << endl;
    }

    return 0;
}