#include <iostream>
#include <sciplot/sciplot.hpp>
#include <vector>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"
#include "../../shared/typedefs.h"

using namespace std;
using namespace capd;
using namespace sciplot;

#define ALG_LOGGER Logger<ProgressIndication>

DTimeMap::SolutionCurve integrateSolution(DMap &map, DVector &point, int order, double initTime, double finalTime)
{
    DOdeSolver solver(map, order);
    DTimeMap timeMap(solver);

    DTimeMap::SolutionCurve solution(initTime);
    timeMap(finalTime, point, solution);
    return solution;
}

Plot2D createPlot()
{
    Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.palette("set1");

    double range = 1.2;
    plot.xrange(-range, range);
    plot.yrange(-range, range);
    plot.size(400, 400);

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

void addPointToPlot(Plot2D &plot, const DVector & point, string label)
{
    plot.drawDots(vector<double>({point[0]}), vector<double>({point[1]})).label(label);
}

void addPointToPlot(Plot2D &plot, const CVector & point, string label)
{
    addPointToPlot(plot, DVector({point[0].real(), point[1].real()}), label);
}

void addCurveToPlot(Plot2D &plot, const vector<double> &x, const vector<double> &y, string label)
{
    plot.drawCurve(x, y).label(label);
}

void showAndSavePlot(Plot2D &plot)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/periodic-orbit/plot.pdf");
}

bool isNearL4(const DVector &point, double epsilon)
{
    DVector L4({0, 0.866025403784438646763723170753});
    return (abs(point[0] - L4[0]) < epsilon && abs(point[1] - L4[1]) < epsilon);
}

int main()
{
    DVector initialPoint({0, -0.576003306698015024546, -0.451540030090179645464, 0 });
    CVector initialPointC({initialPoint[0], initialPoint[1], initialPoint[2], initialPoint[3]});

    double epsilon = 0.01; // normal form will be used on initialPoint + [-epsilon, epsilon]x[-epsilon, epsilon]

    int methodDegree = 20;
    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();
    cout << "normal form imported" << endl;

    DVector point(initialPoint);
    double finalTime = 31; // max time in integration
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, point, 10, 0, finalTime);

    int N = 1000; // numer of samples in a unit of time
    double dt = 1. / N;

    auto plot = createPlot();
    addPointToPlot(plot, initialPoint, "initial point");

    vector<double> solverX, solverY, solverX2, solverY2;
    solverX.reserve(N);
    solverY.reserve(N);
    solverX2.reserve(N);
    solverY2.reserve(N);

    double t = 0;
    CVector lastPoint;
    for( ; t <= finalTime; t += dt)
    {
        if(isNearL4(intSolution(t), epsilon))
        {
            lastPoint = CVector({intSolution(t)[0], intSolution(t)[1], intSolution(t)[2], intSolution(t)[3]});
            break;
        }

        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        solverX2.push_back(-intSolution(t)[0]); // symmetry with respect to mass
        solverY2.push_back(intSolution(t)[1]); 
    }
    cout << "t: " << t << endl;
    addCurveToPlot(plot, solverX, solverY, "integrated solution");
    addCurveToPlot(plot, solverX2, solverY2, "symmetric trajectory");

    vector<double> nfX, nfY;
    for( t = 0; ; t += dt)
    {
        CVector normalFormSolC = normalForm.solution(t, testCase.diagonalization.toDiag(lastPoint));
        normalFormSolC = testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormSolC));
        DVector normalFormSol({normalFormSolC[0].real(), normalFormSolC[1].real()});

        if(!isNearL4(normalFormSol, epsilon))
            break;

        nfX.push_back(normalFormSol[0]);
        nfY.push_back(normalFormSol[1]);
    }
    cout << "t: " << t << endl;
    cout << "diffLeft:\n " <<  (nfX.front() - solverX.back()) << " " << (nfY.front() - solverY.back()) << endl;
    cout << "diffRight:\n" << (nfX.back() - solverX2.back()) << " " << (nfY.back() - solverY2.back()) << endl;
    addCurveToPlot(plot, nfX, nfY, "normal form solution");

    showAndSavePlot(plot);

    return 0;
}