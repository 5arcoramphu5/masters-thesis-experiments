#include <iostream>
#include <fstream>
#include <sciplot/sciplot.hpp>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

#define INIT_TIME 0
#define FINAL_TIME 5

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &point, int order)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(INIT_TIME);
    timeMap(FINAL_TIME, point, solution);
    return solution;
}

Plot2D initializePlot()
{
    Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.palette("set1");

    double range = 2;
    plot.xrange(-range, range);
    plot.yrange(0.866025403784438646763723170753-range, 0.866025403784438646763723170753+range);
    plot.size(1000, 1000);

    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // drawing libration points
    vector<double> librationPointsX({0, 0, 0, 1.19841, -1.19841});
    vector<double> librationPointsY({0, 0.866025403784438646763723170753, -0.866025403784438646763723170753, 0, 0});
    plot.drawDots(librationPointsX, librationPointsY).label("libration points");

    return plot;
}

void addCurveToPlot(Plot2D &plot, const vector<double> &x, const vector<double> &y, string label)
{
    plot.drawCurve(x, y).label(label);
}

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/comparing-trajectories-test/" + filename);
}

// phi = id - Q
// phi^{-1} = id + Q + Q^2 + Q^3 + ... + Q^deg
CVector inverse(const Polynomial<Complex> &phi, const CVector &x, int deg)
{
    CVector result = x;
    CVector curr = x;
    for(int i = 0; i < deg; ++i)
    {
        curr = curr - phi(curr); // Q(curr)
        result += curr;
    }
    return result;
}

void performTest(LDVector &realOriginalPoint, const CVector &normalFormPoint, const PseudoNormalForm &normalForm, TestCasesCollection &testCases)
{
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, realOriginalPoint, 30);

    auto plot = initializePlot();

    int N = 10000;
    double dt = double(FINAL_TIME - INIT_TIME) / N;

    vector<double> intX, intY, normalFormX, normalFormY;
    intX.reserve(N); 
    intY.reserve(N); 
    normalFormX.reserve(N);
    normalFormY.reserve(N);

    for(double t = INIT_TIME; t <= FINAL_TIME; t += dt)
    {
        auto normalFormSol = normalForm.solution(t, normalFormPoint);
        normalFormSol = testCases.PCR3BP_L4.diagonalization.toOriginal(normalForm.getPhi()(normalFormSol));

        cout << "t: " << t << " " << "diff: ";
        for(int j = 0; j < 4; ++j)
            cout << ((Complex)intSolution(t)[j] - normalFormSol[j]) << " ";
        cout << endl;

        intX.push_back(intSolution(t)[0]);
        intY.push_back(intSolution(t)[1]);
        normalFormX.push_back(normalFormSol[0].real());
        normalFormY.push_back(normalFormSol[1].real());
    }

    addCurveToPlot(plot, intX, intY, "integrated solution");
    addCurveToPlot(plot, normalFormX, normalFormY, "normal form solution");
    showAndSavePlot(plot, "plot.pdf");
}

int main()
{
    int methodDegree = 15;
    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();

    // Complex lambda1(0.0000001, 0), lambda2(-0.0000001, 0); // the closer to 0 the closer to p transformedPoint will be
    // CVector normalFormPoint({ lambda1, lambda2, conj(lambda1), conj(lambda2) });
    // auto originalSpacePoint =  testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormPoint));

    // DVector realOriginalPoint(4);
    // for(int i = 0; i < 4; ++i)  realOriginalPoint[i] = originalSpacePoint[i].real();

    // performTest(realOriginalPoint, normalFormPoint, normalForm, testCases);

    LDVector realOriginalPoint({0, 0.9, 0.1, 0});
    CVector originalPoint ({realOriginalPoint[0], realOriginalPoint[1], realOriginalPoint[2], realOriginalPoint[3]});
    CVector normalFormPoint = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(originalPoint), 15);
    performTest(realOriginalPoint, normalFormPoint, normalForm, testCases);

    return 0;
}