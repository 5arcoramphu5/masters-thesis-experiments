#include <iostream>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"
#include "../../shared/typedefs.h"

using namespace std;
using namespace capd;

#define ALG_LOGGER Logger<None, MathematicaFormatPolynomialPrinting, 5>

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

void printSolution(DTimeMap::SolutionCurve &s, int N)
{
    double dt = double(FINAL_TIME - INIT_TIME) / N;
    for(double t = INIT_TIME; t <= FINAL_TIME; t += dt)
    {
        cout << "t: " << t << ",\t solution: " << s(t) << endl;
    }
}

int main()
{
    TestCasesCollection testCases(10);

    NormalFormFinder<ALG_LOGGER> finder(4,  testCases.PCR3BP);
    auto normalForm = finder.calculatePseudoNormalForm();

    Complex lambda1(0.000001, 0), lambda2(-0.000001, 0); // the closer to 0 the closer to p transformedPoint will be
    CVector point({ lambda1, lambda2, conj(lambda1), conj(lambda2) });
    auto transformedPoint =  testCases.PCR3BP.toOriginal(normalForm.getPhi()(point));
    cout << "transformed point:" << transformedPoint << "\n\n";

    DVector realTransformedPoint(4);
    for(int i = 0; i < 4; ++i)
        realTransformedPoint[i] = transformedPoint[i].real();

    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, realTransformedPoint, 20);
    printSolution(intSolution, 10);

    return 0;
}