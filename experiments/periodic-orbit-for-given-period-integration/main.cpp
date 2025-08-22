#include <iostream>
#include <iomanip>
#include <vector>
#include <sciplot/sciplot.hpp>

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/lib.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace capd::matrixAlgorithms;
using namespace sciplot;

#define PLOT_SCALE 1000
#define L4_Y 0.866025403784438646763723170753
#define PLOT_SIZE_PIXELS 300

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

    plot.xrange(-1, 1);
    plot.yrange(-0.7, 1.3);
    plot.size(PLOT_SIZE_PIXELS, PLOT_SIZE_PIXELS);

    double scale = 10*PLOT_SCALE;
    plotCloseUp.xrange(-epsilon*scale, epsilon*scale);
    plotCloseUp.yrange(0.866025403784438646763723170753-epsilon*scale, 0.866025403784438646763723170753+epsilon*scale);
    plotCloseUp.size(PLOT_SIZE_PIXELS, PLOT_SIZE_PIXELS);

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
    plot.drawPoints(librationPointsX, librationPointsY).label("libration points");
    plotCloseUp.drawPoints( vector<double>({librationPointsX[1]}), vector<double>({librationPointsY[1]})).label("L4");
}

void addCurveToPlots(Plot2D &plot, Plot2D &plotCloseUp, const vector<double> &x, const vector<double> &y, string label)
{
    plot.drawCurve(x, y).label(label);

    int n = x.size();
    vector<long double> x_scaled(n), y_scaled(n);
    for(int i = 0; i < n; ++i)
    {
        x_scaled[i] = PLOT_SCALE * x[i];
        y_scaled[i] = PLOT_SCALE * (y[i] - L4_Y) + L4_Y;
    }

    plotCloseUp.drawCurve(x_scaled, y_scaled).label(label);
}

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(PLOT_SIZE_PIXELS, PLOT_SIZE_PIXELS);

    canvas.show();
    canvas.save("experiments/periodic-orbit-for-given-period-integration/" + filename);
}

void showAndSavePlots(Plot2D &plot, Plot2D &plotCloseUp)
{
    showAndSavePlot(plot, "plot.pdf");
    showAndSavePlot(plotCloseUp, "plotCloseUp.pdf");
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

// otrzymuje (y, v_x) 
// oblicza punkt końcowy dla casu halfPeriod przy punkcie początkowym (0, y, v_x, 0)
// zwraca rzutowanie punktu końcowego na zmienne (x, v_y)
// celem jest znalezienie miejsca zerowego, więc punktu (0, y', v_x', 0)
// zwraca wynik i pochodną
pair<LDVector, LDMatrix> computeF(const LDVector &initialPoint, double halfPeriod, LDMap &map, int order)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDVector point({0, initialPoint[0], initialPoint[1], 0});
    LDMatrix der;
    timeMap(halfPeriod, point, der);

    LDVector result({point[0], point[3]});
    
    return make_pair(result, der);
}

LDVector getInitialPoint(double halfPeriod, const LDVector &approximation, LDMap &map, double eps)
{
    int solverOrder = 15;
    LDVector point(approximation);

    while(true)
    {
        auto [f, f_der] = computeF(point, halfPeriod, map, solverOrder);

        if(abs(f[0]) < eps && abs(f[1]) < eps)
            break;

        cout << "f: " << f << endl;

        // trunceting matrix to take projections into account
        LDMatrix der({ {f_der[0][1], f_der[0][2]}, {f_der[3][1], f_der[3][2]}});

        // inverse of 2x2 matrix
        long double det = der[0][0]*der[1][1] - der[0][1]*der[1][0];
        LDMatrix invDer({   {der[1][1]/det, -der[0][1]/det}, 
                            {-der[1][0]/det, der[0][0]/det} });
        // invDer * f
        LDVector invDerF({invDer[0][0]*f[0] + invDer[0][1]*f[1], invDer[1][0]*f[0] + invDer[1][1]*f[1] });
        point =  point - invDerF;
    }

    return LDVector({0, point[0], point[1], 0});
}

int main(int argc, char* argv[])
{
    if(argc != 3) throw runtime_error("wrong number of arguments");

    double period = atof(argv[1]);
    double eps = atof(argv[2]);

    cout << std::setprecision(19);

    TestCasesCollection testCases(20);
    auto &testCase = testCases.PCR3BP_L4;

    LDVector approximation({-0.576003306698015024546, -0.451540030090179645464}); // w przybliżeniu orbita okresowa zaczyna się w (0, x, y, 0)
    auto initialPoint = getInitialPoint(period/2, approximation, testCases.PCR3BP_dmap, eps);
    cout << "initial point: " << initialPoint << endl;

    double finalTime = period / 2; // max time in integration
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, initialPoint, 10, 0, finalTime);

    int N = 1000; // numer of samples in a unit of time
    double intDt = 1. / N;

    vector<double> solverX, solverY, solverX2, solverY2;
    double epsilon = 1; 

    for(double t = 0; t <= finalTime; t += intDt)
    {
        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        // symmetry with respect to mass
        solverX2.push_back(-intSolution(t)[0]); 
        solverY2.push_back(intSolution(t)[1]); 

        LDVector diff({intSolution(t)[0], intSolution(t)[1] - 0.866025403784438646763723170753});
        if(diff.euclNorm() < epsilon)
            epsilon = diff.euclNorm();
    }
    
    Plot2D plot, plotCloseUp;
    initializePlots(plot, plotCloseUp, epsilon);
    addCurveToPlots(plot, plotCloseUp, solverX, solverY, "integrated solution");
    addCurveToPlots(plot, plotCloseUp, solverX2, solverY2, "symmetric trajectory");
    showAndSavePlots(plot, plotCloseUp);

    return 0;
}