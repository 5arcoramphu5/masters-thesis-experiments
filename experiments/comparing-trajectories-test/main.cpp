#include <iostream>
#include <vector>
#include <fstream>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"
#include "../../shared/typedefs.h"

using namespace std;
using namespace capd;

#define ALG_LOGGER Logger<None>
#define LOGGER Logger<Minimal, MathematicaFormatPolynomialPrinting, 10>

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

void compareSolutions(DTimeMap::SolutionCurve &s1, const vector<CVector> &s2, int N)
{
    double dt = double(FINAL_TIME - INIT_TIME) / N;
    int i = 0;
    for(double t = INIT_TIME; t <= FINAL_TIME; t += dt)
    {
        cout << "t: " << t << ",\t s1: " << s1(t) << ",\t s2: " << s2[i] << endl;
        cout << "diff: ";
            for(int j = 0; j < 4; ++j)
                cout << (s1(t)[j] - s2[i][j].real()) << " ";
        cout << endl;
        i++;
    }
}

int main()
{
    int methodDegree = 7;
    TestCasesCollection testCases(10);
    string path = "experiments/comparing-trajectories-test/results/";

    ifstream fileN(path + "N_deg" + to_string(methodDegree) + ".txt");
    auto N = Polynomial<Complex>::deserialize(fileN);
    fileN.close();
    LOGGER::log<Minimal>("N:\n", N);

    ifstream filePhi(path + "Phi_deg" + to_string(methodDegree) + ".txt");
    auto phi = Polynomial<Complex>::deserialize(filePhi);
    filePhi.close();

    Complex lambda1(0.000001, 0), lambda2(-0.000001, 0); // the closer to 0 the closer to p transformedPoint will be
    CVector point({ lambda1, lambda2, conj(lambda1), conj(lambda2) });
    auto transformedPoint =  testCases.PCR3BP.toOriginal(phi(point));
    cout << "transformed point: " << transformedPoint << "\n\n";

    DVector realTransformedPoint(4);
    for(int i = 0; i < 4; ++i)
        realTransformedPoint[i] = transformedPoint[i].real();

    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, realTransformedPoint, 20);

    // from mathematica
    vector<CVector> normal_form_trajectory = {CVector({Complex(0.000001, 0.), Complex(-0.000001, 0.), Complex(0.000001, 0.), Complex(-0.000001, 0.)}), CVector({Complex(0.00000111412, 0.000000213964), Complex(-0.000000865645, 0.000000166245), Complex(0.00000111412, -0.000000213964), Complex(-0.000000865645, -0.000000166245)}), CVector({Complex(0.0000011957, 0.000000476789), Complex(-0.000000721569, 0.000000287752), Complex(0.0000011957, -0.000000476789), Complex(-0.000000721569, -0.000000287752)}), CVector({Complex(0.00000123047, 0.000000787143), Complex(-0.000000576654, 0.000000368929), Complex(0.00000123047, -0.000000787143), Complex(-0.000000576654, -0.000000368929)}), CVector({Complex(0.00000120282, 0.00000114044), Complex(-0.000000437751, 0.000000415104), Complex(0.00000120282, -0.00000114044), Complex(-0.000000437751, -0.000000415104)}), CVector({Complex(0.00000109641, 0.00000152823), Complex(-0.000000309867, 0.000000431984), Complex(0.00000109641, -0.00000152823), Complex(-0.000000309867, -0.000000431984)})};

    for(int i = 0; i < normal_form_trajectory.size(); ++i)
        normal_form_trajectory[i] = testCases.PCR3BP.toOriginal(phi(normal_form_trajectory[i]));

    compareSolutions(intSolution, normal_form_trajectory, 5);

    return 0;
}