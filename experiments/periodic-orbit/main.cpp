#include <iostream>
#include <sciplot/sciplot.hpp>
#include <vector>

#include "capd/capdlib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

#define ALG_LOGGER Logger<ProgressIndication>

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &initialPoint, int order, double initTime, double finalTime)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(initTime);
    LDVector point(initialPoint);
    timeMap(finalTime, point, solution);
    return solution;
}

void initializePlots(Plot2D &plot, Plot2D &plotCloseUp, double epsilon)
{
    plot.xlabel("x");
    plot.ylabel("y");
    plotCloseUp.xlabel("x");
    plotCloseUp.ylabel("y");

    plot.palette("set1");
    plotCloseUp.palette("set1");

    double range = 1.2;
    plot.xrange(-range, range);
    plot.yrange(-range, range);
    plot.size(1000, 1000);

    double scale = 3;
    plotCloseUp.xrange(-epsilon*scale, epsilon*scale);
    plotCloseUp.yrange(0.866025403784438646763723170753-epsilon*scale, 0.866025403784438646763723170753+epsilon*scale);
    plotCloseUp.size(1000, 1000);

    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);
    plotCloseUp.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // drawing libration points
    vector<double> librationPointsX({0, 0, 0, 1.19841, -1.19841});
    vector<double> librationPointsY({0, 0.866025403784438646763723170753, -0.866025403784438646763723170753, 0, 0});
    plot.drawDots(librationPointsX, librationPointsY).label("libration points");
    plotCloseUp.drawDots( vector<double>({librationPointsX[1]}), vector<double>({librationPointsY[1]})).label("L4");
}

void addPointToPlots(Plot2D &plot, Plot2D &plotCloseUp, const DVector & point, string label)
{
    plot.drawDots(vector<double>({point[0]}), vector<double>({point[1]})).label(label);
    plotCloseUp.drawDots(vector<double>({point[0]}), vector<double>({point[1]})).label(label);
}

void addPointToPlots(Plot2D &plot, Plot2D &plotCloseUp, const CVector & point, string label)
{
    addPointToPlots(plot, plotCloseUp, DVector({point[0].real(), point[1].real()}), label);
}

void addCurveToPlots(Plot2D &plot, Plot2D &plotCloseUp, const vector<double> &x, const vector<double> &y, string label)
{
    plot.drawCurve(x, y).label(label);
    plotCloseUp.drawCurve(x, y).label(label);
}

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/periodic-orbit/" + filename);
}

void showAndSavePlots(Plot2D &plot, Plot2D &plotCloseUp)
{
    showAndSavePlot(plot, "plot.pdf");
    showAndSavePlot(plotCloseUp, "plotCloseUp.pdf");
}

bool isNearL4(const LDVector &point, double epsilon)
{
    DVector L4({0, 0.866025403784438646763723170753});
    return (abs(point[0] - L4[0]) < epsilon && abs(point[1] - L4[1]) < epsilon);
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

int main()
{
    double period = 30.4301813921308539346;
    LDVector initialPoint({0, -0.5760033133504633678, -0.4515400305608129141, 0});
    CVector initialPointC({initialPoint[0], initialPoint[1], initialPoint[2], initialPoint[3]});

    double epsilon = 0.001; // normal form will be used on initialPoint + [-epsilon, epsilon]x[-epsilon, epsilon]

    int methodDegree = 20;
    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();
    cout << "normal form imported" << endl;

    double finalTime = period / 2; // max time in integration
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, initialPoint, 10, 0, finalTime);

    int N = 1000; // numer of samples in a unit of time
    double intDt = 1. / N;

    Plot2D plot, plotCloseUp;
    initializePlots(plot, plotCloseUp, epsilon);

    vector<double> solverX, solverY, solverX2, solverY2;
    solverX.reserve(N);
    solverY.reserve(N);
    solverX2.reserve(N);
    solverY2.reserve(N);

    CVector lastPoint;
    double timeLeft;
    for(double t = 0; t <= finalTime; t += intDt)
    {
        if(isNearL4(intSolution(t), epsilon))
        {
            lastPoint = CVector({intSolution(t)[0], intSolution(t)[1], intSolution(t)[2], intSolution(t)[3]});
            timeLeft = finalTime-t;
            break;
        }

        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        // symmetry with respect to mass
        solverX2.push_back(-intSolution(t)[0]); 
        solverY2.push_back(intSolution(t)[1]); 
    }

    addCurveToPlots(plot, plotCloseUp, solverX, solverY, "integrated solution");
    addCurveToPlots(plot, plotCloseUp, solverX2, solverY2, "symmetric trajectory");

    vector<double> nfX, nfY;
    double nfDt = 0.1;
    auto newInitialPoint = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(lastPoint), 100);

    for(double t = 0; t < 2*timeLeft; t += nfDt)
    {
        CVector normalFormSolC = normalForm.solution(t, newInitialPoint);
        normalFormSolC = testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormSolC));
        DVector normalFormSol({normalFormSolC[0].real(), normalFormSolC[1].real()});

        nfX.push_back(normalFormSol[0]);
        nfY.push_back(normalFormSol[1]);
    }

    cout << "diffLeft:\n " <<  (nfX.front() - solverX.back()) << " " << (nfY.front() - solverY.back()) << endl;
    cout << "diffRight:\n" << (nfX.back() - solverX2.back()) << " " << (nfY.back() - solverY2.back()) << endl;

    addCurveToPlots(plot, plotCloseUp, nfX, nfY, "normal form solution");
    showAndSavePlots(plot, plotCloseUp);

    return 0;
}