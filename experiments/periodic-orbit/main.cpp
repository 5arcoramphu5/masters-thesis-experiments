#include <iostream>
#include <sciplot/sciplot.hpp>
#include <vector>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &initialPoint, int order, double initTime, double finalTime)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(initTime);
    LDVector point(initialPoint);
    timeMap(finalTime, point, solution);
    return solution;
}

Plot2D initializePlot()
{
    Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");

    plot.palette("set1");

    double range = 1.2;
    plot.xrange(-range, range);
    plot.yrange(-range, range);
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

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/periodic-orbit/" + filename);
}

int main()
{
    double period = 30.4301813921308539346*2;
    LDVector initialPoint({0, -0.576003306698015024546, -0.451540030090179645464, 0});
    CVector initialPointC({initialPoint[0], initialPoint[1], initialPoint[2], initialPoint[3]});

    TestCasesCollection testCases(20);
    auto &testCase = testCases.PCR3BP_L4;

    double finalTime = period / 2; // max time in integration
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, initialPoint, 15, 0, finalTime);

    cout << intSolution(finalTime) << endl;

    int N = 100; // numer of samples in a unit of time
    double intDt = 1. / N;

    auto plot = initializePlot();
    vector<double> solverX, solverY, solverX2, solverY2;

    for(double t = 0; t <= finalTime; t += intDt)
    {
        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        // symmetry with respect to mass
        solverX2.push_back(-intSolution(t)[0]); 
        solverY2.push_back(intSolution(t)[1]); 
    }

    plot.drawCurve(solverX, solverY).label("integrated solution");
    plot.drawCurve(solverX2, solverY2).label("symmetric trajectory");
    showAndSavePlot(plot, "plot.pdf");

    return 0;
}