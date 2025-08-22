#include <iostream>
#include <iomanip>
#include <vector>
#include <sciplot/sciplot.hpp>

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/lib.h"

#include "../../shared/calculations/calculations.h"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace capd::matrixAlgorithms;
using namespace sciplot;

#define PLOT_SCALE 1000
#define L4_Y 0.866025403784438646763723170753
#define SOLVER_ORDER 10

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, const LDVector &initialPoint, double initTime, double finalTime)
{
    LDOdeSolver solver(map, SOLVER_ORDER);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(initTime);
    LDVector point(initialPoint);
    timeMap(finalTime, point, solution);
    return solution;
}

void initializePlots(Plot2D &plot, Plot2D &plotCloseUp, LDVector lastPointInt)
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

    LDVector L4({0, 0.866025403784438646763723170753});
    LDVector lastPointInt2D({lastPointInt[0], lastPointInt[1]});
    double epsilon = (lastPointInt2D - L4).euclNorm();
    double rangeMul = 10*PLOT_SCALE;
    plotCloseUp.xrange(-epsilon*rangeMul, epsilon*rangeMul);
    plotCloseUp.yrange(L4_Y-epsilon*rangeMul, L4_Y+epsilon*rangeMul);
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

void showAndSavePlot(Plot2D &plot, string filename, string plotTitle)
{
    Figure fig({{plot}});
    fig.title(plotTitle);
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/periodic-orbit-for-given-period-normal-form-2/" + filename);
}

void showAndSavePlots(Plot2D &plot, Plot2D &plotCloseUp, string plotTitle)
{
    showAndSavePlot(plot, "plot.pdf", plotTitle);
    showAndSavePlot(plotCloseUp, "plotCloseUp.pdf", plotTitle);
}

LDVector normalFormSolution(double t, const LDVector &point, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    CVector pointC({point[0], point[1], point[2], point[3]});
    auto nfCoordsPoint = inverse(normalForm.getPhi(), diag.toDiag(pointC));

    CVector normalFormSol = normalForm.solution(t, nfCoordsPoint);
    normalFormSol = diag.toOriginal(normalForm.getPhi()(normalFormSol));

    return LDVector({normalFormSol[0].real(), normalFormSol[1].real(), normalFormSol[2].real(), normalFormSol[3].real()});
}

pair<LDVector, LDMatrix> integrationSolution(double t, const LDVector &initialPoint, LDMap &map)
{
    LDOdeSolver solver(map, SOLVER_ORDER);
    LDTimeMap timeMap(solver);

    LDVector point(initialPoint);
    LDMatrix intDer;
    point = timeMap(t, point, intDer);

    return make_pair(point, intDer);
}

// otrzymuje (y, v_x, x*, y*, v_x*, v_y*)
// 
// celem jest znalezienie miejsca zerowego, więc punktu (0, y', v_x', 0)
// zwraca wynik i pochodną
pair<LDVector, LDMatrix> computeF(const LDVector &point, double halfPeriod, double intTime, LDMap &map, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    // compute value of the function
    LDVector result(6);

    LDVector nfSol = normalFormSolution(halfPeriod - intTime, LDVector({point[2], point[3], point[4], point[5]}), normalForm, diag);
    auto [intSol, intDer] = integrationSolution(intTime, LDVector({0, point[0], point[1], 0}), map);

    // cout << "nfSol: " << nfSol << endl;
    // cout << "intSol: " << intSol << endl;

    for(int i = 0; i < 4; ++i)
        result[i] = point[i+2] - intSol[i];

    result[4] = nfSol[0];
    result[5] = nfSol[3];

    // compute derivative
    LDMatrix der(6, 6);
    auto nfDer = normalFormOriginalSolutionDerivative(halfPeriod - intTime, LDVector({point[2], point[3], point[4], point[5]}), normalForm, diag);

    der[4][0] = der[4][1] = der[5][0] = der[5][1] = 0;

    for(int i = 0; i < 4; ++i)
    {
        der[i][2+i] = 1;

        der[i][0] = -intDer[i][1]; 
        der[i][1] = -intDer[i][2];
        
        der[4][2+i] = nfDer[0][i];
        der[5][2+i] = nfDer[3][i];
    }
    
    return make_pair(result, der);
}

pair<LDVector, LDVector> getNormalFormInitialPoint(double halfPeriod, double intTime, const LDVector initialPoint, LDMap &map, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag, double eps)
{
    auto afterInt = integrationSolution(intTime, initialPoint, map).first;

    // (y, dx, x1, y1, dx1, dx2)
    LDVector point({initialPoint[1], initialPoint[2], afterInt[0], afterInt[1], afterInt[2], afterInt[3]});
    int maxIters = 20;

    auto bestPoint = point;
    optional<LDVector> bestF = std::nullopt;

    for(int i = 0; i < maxIters; ++i)
    {
        cout << "i: " << i << endl;
        auto [f, f_der] = computeF(point, halfPeriod, intTime, map, normalForm, diag);

        if(!bestF.has_value() || abs(f) < abs(bestF.value()))
        {
            bestF = f;
            bestPoint = point;
        }
        
        if(abs(f[0]) < eps && abs(f[1]) < eps && abs(f[2]) < eps && abs(f[3]) < eps && abs(f[4]) < eps && abs(f[5]) < eps)
        {
            cout << "found initial point: " << point << endl;
            cout << "value: " << f << endl; 
            break;
        }
        cout << f << endl;
        
        // invDerF = (f_der)^{-1} * f
        // f_der * invDerF = f - rownanie liniowe
        LDVector invDerF = gauss(f_der, f);

        point = point - invDerF;
    }

    return make_pair(LDVector({0, bestPoint[0], bestPoint[1], 0}), LDVector({bestPoint[2], bestPoint[3], bestPoint[4], bestPoint[5]}));
}

int main(int argc, char* argv[])
{
    if(argc != 5) throw runtime_error("wrong number of arguments");

    double period = atof(argv[1]);
    double intTime = atof(argv[2]);
    double eps = atof(argv[3]);
    int methodDegree = atof(argv[4]);

    if(intTime >= period/2)
        throw runtime_error("integration time must be less than half of the period");

    double timeLeft = period/2 - intTime;

    cout << std::setprecision(19);

    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();
    cout << "normal form imported" << endl;

    LDVector approxIntInitialPoint({0, -0.576003306698015024546, -0.451540030090179645464, 0});

    auto[intInitialPoint, nfInitialPoint] = getNormalFormInitialPoint(period/2, intTime, approxIntInitialPoint, testCases.PCR3BP_dmap, normalForm, testCase.diagonalization, eps);
    cout << "intInitialPoint: " << intInitialPoint << endl;
    cout << "nfInitialPoint: " << nfInitialPoint << endl;

    // plot integrated solution
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, intInitialPoint, 0, intTime);

    int N = 100; // numer of samples in a unit of time
    double intDt = 1. / N;

    vector<double> solverX, solverY, solverX2, solverY2;
    for(double t = 0; t <= intTime; t += intDt)
    {
        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        // symmetry with respect to mass
        solverX2.push_back(-intSolution(t)[0]);
        solverY2.push_back(intSolution(t)[1]);
    }

    // plot normal form solution
    vector<double> nfX, nfY;
    double nfDt = timeLeft / N;

    CVector nfInitialPointC({nfInitialPoint[0], nfInitialPoint[1], nfInitialPoint[2], nfInitialPoint[3]});
    auto nfInitialPointDiag = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(nfInitialPointC));

    for(double t = 0; t <= timeLeft; t += nfDt)
    {
        CVector nfSol = normalForm.solution(t, nfInitialPointDiag);
        CVector nfSolOrigCoords = testCase.diagonalization.toOriginal(normalForm.getPhi()(nfSol));

        nfX.push_back(nfSolOrigCoords[0].real());
        nfY.push_back(nfSolOrigCoords[1].real());
    }

    cout << "diffLeft:\n " <<  (nfX.front() - solverX.back()) << " " << (nfY.front() - solverY.back()) << endl;
    cout << "diffRight:\n" << (nfX.back() - solverX2.back()) << " " << (nfY.back() - solverY2.back()) << endl;

    Plot2D plot, plotCloseUp;
    initializePlots(plot, plotCloseUp, intSolution(intTime));

    addCurveToPlots(plot, plotCloseUp, solverX, solverY, "integrated solution");
    addCurveToPlots(plot, plotCloseUp, solverX2, solverY2, "symmetric trajectory");
    addCurveToPlots(plot, plotCloseUp, nfX, nfY, "normal form solution");

    string title = "period: " + to_string(period) + ", int time: " + to_string(intTime) + ", epsilon: " + to_string(eps);
    showAndSavePlots(plot, plotCloseUp, title);

    return 0;
}