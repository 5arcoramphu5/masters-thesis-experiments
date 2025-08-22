#include <iostream>
#include <fstream>
#include <sciplot/sciplot.hpp>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/calculations/calculations.h"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

#define INIT_TIME 0
#define FINAL_TIME 4

#define L4_Y 0.866025403784438646763723170753 // y coordinate of L4 point in PCR3BP

#define PLOT_SIZE_PIXELS 300

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &point, int order)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(INIT_TIME);
    timeMap(FINAL_TIME, point, solution);
    return solution;
}

Plot2D initializePlot(double range)
{
    Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.palette("set1");

    plot.xrange(-range, range);
    plot.yrange(L4_Y-range, L4_Y+range);
    plot.size(PLOT_SIZE_PIXELS, PLOT_SIZE_PIXELS);

    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // drawing libration points
    vector<double> librationPointsX({0, 0, 0, 1.19841, -1.19841});
    vector<double> librationPointsY({0, L4_Y, -L4_Y, 0, 0});
    plot.drawPoints(librationPointsX, librationPointsY).label("libration points");

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
    canvas.size(PLOT_SIZE_PIXELS, PLOT_SIZE_PIXELS);

    canvas.show();
    canvas.save("experiments/comparing-trajectories-test/" + filename);
}

void performTest(LDVector &realOriginalPoint, const CVector &normalFormPoint, const PseudoNormalForm &normalForm, TestCasesCollection &testCases, double plotRange)
{
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, realOriginalPoint, 30);

    auto plot = initializePlot(plotRange);

    int N = 10000;
    double dt = double(FINAL_TIME - INIT_TIME) / N;

    double epsilon = 0.01; // difference between solutions in shown for points with distance from L4 <= epsilon

    vector<double> intX, intY, normalFormX, normalFormY;

    double maxDiff = 0;
    for(double t = INIT_TIME; t <= FINAL_TIME; t += dt)
    {
        auto normalFormSol = normalForm.solution(t, normalFormPoint);
        normalFormSol = testCases.PCR3BP_L4.diagonalization.toOriginal(normalForm.getPhi()(normalFormSol));

        if(abs(intSolution(t)[0]) <= epsilon && abs(intSolution(t)[1] - L4_Y) <= epsilon)
        {
            cout << "t: " << t << " " << "diff: ";
            for(int j = 0; j < 4; ++j)
            {
                auto diff = ((Complex)intSolution(t)[j] - normalFormSol[j]);
                cout << diff << " ";

                if(abs(diff) > maxDiff)
                    maxDiff = abs(diff);
            }
            cout << endl;
        }

        intX.push_back(intSolution(t)[0]);
        intY.push_back(intSolution(t)[1]);
        normalFormX.push_back(normalFormSol[0].real());
        normalFormY.push_back(normalFormSol[1].real());
    }
    cout << "max diff: " << maxDiff << endl;

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

    // LDVector realOriginalPoint({0, 0.9, 0.1, 0});
    // double plotRange = 1;

    // LDVector realOriginalPoint({0.001, 0.865, 0.001, 0});
    // double plotRange = 0.05;

    LDVector realOriginalPoint({0.00001, 0.86601, 0.00001, 0});
    double plotRange = 0.0005;
    
    CVector originalPoint ({realOriginalPoint[0], realOriginalPoint[1], realOriginalPoint[2], realOriginalPoint[3]});
    CVector normalFormPoint = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(originalPoint));
    performTest(realOriginalPoint, normalFormPoint, normalForm, testCases, plotRange);

    return 0;
}