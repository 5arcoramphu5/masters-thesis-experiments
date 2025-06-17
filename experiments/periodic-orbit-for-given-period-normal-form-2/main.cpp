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

void addCurveToPlots(Plot2D &plot, Plot2D &plotCloseUp, const vector<double> &x, const vector<double> &y, string label)
{
    plot.drawCurve(x, y).label(label);
    plotCloseUp.drawCurve(x, y).label(label);
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

CMatrix inverseDer(const Polynomial<Complex> &phi, const CVector &x, int deg)
{
    CMatrix result = CMatrix::Identity(4);
    CMatrix currDer = CMatrix::Identity(4);
    CVector currVal = x;

    for(int i = 0; i < deg; ++i)
    {
        currVal = currVal - phi(currVal); // Q(curr)
        currDer = (-phi.derivative(currVal) + CMatrix::Identity(4)) * currDer; // Df(Q(curr)) * Df(f^{-1}(x))
        result += currDer;
    }
    return result;
}


CMatrix normalFormSolutionDerivative(double t, CVector point, const PseudoNormalForm &normalForm)
{
    CVector mulV({ point[0]*point[1], point[2]*point[3]});
    CMatrix derMul({ // derivative of (x[0]*x[1], x[2]*x[3])
        { point[1], point[0], 0, 0},
        { 0, 0, point[2], point[3]}
    });

    auto a1 = normalForm.geta1()(mulV)[0];
    auto a2 = normalForm.geta2()(mulV)[0]; 

    auto derA1 = normalForm.geta1().derivative(mulV);
    auto derA2 = normalForm.geta2().derivative(mulV);

    auto derA1mul = (derA1 * derMul)[0];
    auto derA2mul = (derA2 * derMul)[0];
    
    auto exp1 = exp(a1*t);
    auto exp2 = exp(-a1*t);
    auto exp3 = exp(a2*t);
    auto exp4 = exp(-a2*t);

    Complex one = 1;

    CMatrix result({
        { (one + derA1mul[0]*point[0]) * exp1, point[0]*exp1*derA1mul[1], point[0]*exp1*derA1mul[2], point[0]*exp1*derA1mul[3] },
        { -point[1]*exp2*derA1mul[0], (one - derA1mul[1]*point[1]) * exp2, -point[1]*exp2*derA1mul[2], -point[1]*exp2*derA1mul[3] },
        { point[2]*exp3*derA2mul[0], point[2]*exp3*derA2mul[1], (one + derA2mul[2]*point[2]) * exp3, point[2]*exp3*derA2mul[3] },
        { -point[3]*exp4*derA2mul[0], -point[3]*exp4*derA2mul[1], -point[3]*exp4*derA2mul[2], (one - derA2mul[3]*point[3]) * exp4}
    });

    return result;
}

LDMatrix normalFormOriginalSolutionDerivative(double t, const LDVector &point, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    CVector pointC({point[0], point[1], point[2], point[3]});
    CVector diagPoint = diag.toDiag(pointC);

    int inverseDeg = 100;

    CVector invPhiDiag = inverse(normalForm.getPhi(), diagPoint, inverseDeg);
    CVector solInvPhiDIag = normalForm.solution(t, invPhiDiag);

    CMatrix derPhi = normalForm.getPhi().derivative(solInvPhiDIag);
    CMatrix derNFSol = normalFormSolutionDerivative(t, invPhiDiag, normalForm);
    CMatrix derInvPhi = inverseDer(normalForm.getPhi(), diagPoint, inverseDeg);

    CMatrix resultComplex = diag.getinvJ() * derPhi * derNFSol * derInvPhi * diag.getJ();

    LDMatrix result(4, 4);
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            result[i][j] = resultComplex[i][j].real();

    return result;
}

// otrzymuje (y, v_x) 
// oblicza punkt końcowy dla casu halfPeriod przy punkcie początkowym (0, y, v_x, 0)
// zwraca rzutowanie punktu końcowego na zmienne (x, v_y)
// celem jest znalezienie miejsca zerowego, więc punktu (0, y', v_x', 0)
// zwraca wynik i pochodną
pair<LDVector, LDMatrix> computeF(const LDVector &initialPoint, double halfPeriod, double intTime, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    // compute value of the function
    CVector initialPointC({initialPoint[0], initialPoint[1], initialPoint[2], initialPoint[3]});
    auto newInitialPoint = inverse(normalForm.getPhi(), diag.toDiag(initialPointC), 100);

    CVector normalFormSol = normalForm.solution(halfPeriod - intTime, newInitialPoint);
    normalFormSol = diag.toOriginal(normalForm.getPhi()(normalFormSol));
    LDVector result({normalFormSol[0].real(), normalFormSol[3].real()});

    // compute derivative
    auto der = normalFormOriginalSolutionDerivative(halfPeriod - intTime, initialPoint, normalForm, diag);
    
    return make_pair(result, der);
}

LDVector getNormalFormInitialPoint(double halfPeriod, double intTime, const LDVector &approximation, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag, double eps)
{
    LDVector point({approximation[1], approximation[2]});
    int maxIters = 20;

    for(int i = 0; i < maxIters; ++i)
    {
        LDVector fullPoint({approximation[0], point[0], point[1], approximation[3]});
        auto [f, f_der] = computeF(fullPoint, halfPeriod, intTime, normalForm, diag);

        if(abs(f[0]) < eps && abs(f[1]) < eps)
        {
            cout << "found initial point: " << point << endl;
            cout << "value: " << f << endl; 
            break;
        }
        cout << f << endl;
        // cout << f_der << endl;

        // trunceting matrix to take projections into account
        LDMatrix der({ 
            {f_der[0][1], f_der[0][2]}, 
            {f_der[3][1], f_der[3][2]} });
        // cout << "der:\n" << der << endl;

        // inverse of 2x2 matrix
        long double det = der[0][0]*der[1][1] - der[0][1]*der[1][0];
        LDMatrix invDer({   {der[1][1]/det, -der[0][1]/det}, 
                            {-der[1][0]/det, der[0][0]/det} });

        // cout << "invDer:\n" << invDer << endl;
        // invDer * f
        LDVector invDerF({invDer[0][0]*f[0] + invDer[0][1]*f[1], invDer[1][0]*f[0] + invDer[1][1]*f[1] });
        // cout << "invDerF: " << invDerF << endl;
        point =  point - invDerF;

        // if(abs(invDerF[0]) < eps && abs(invDerF[1]) < eps)
        // {
        //     cout << "found initial point: " << point << endl;
        //     cout << "value: " << f << endl; 
        //     break;
        // }
    }

    return LDVector({approximation[0], point[0], point[1], approximation[3]});
}

int main(int argc, char* argv[])
{
    if(argc != 4) throw runtime_error("wrong number of arguments");

    double period, intTime, eps;
    period = atof(argv[1]);
    intTime = atof(argv[2]);
    eps = atof(argv[3]);

    if(intTime >= period/2)
        throw runtime_error("integration time must be less than half of the period");

    int methodDegree = 20;
    double timeLeft = period/2 - intTime;

    cout << std::setprecision(19);

    TestCasesCollection testCases(methodDegree+1);
    auto &testCase = testCases.PCR3BP_L4;

    ifstream file("shared/precalculated-normal-forms/PCR3BP_L4_deg" + to_string(methodDegree) + ".txt");
    auto normalForm = PseudoNormalForm::deserialize(file);
    file.close();
    cout << "normal form imported" << endl;

    LDVector initialPoint({0, -0.576003306698015024546, -0.451540030090179645464, 0});
    auto initialSolution = integrateSolution(testCases.PCR3BP_dmap, initialPoint, 10, 0, intTime);

    LDVector normalFormInitialPoint = getNormalFormInitialPoint(period/2, intTime, initialSolution(intTime), normalForm, testCase.diagonalization, eps);
    cout << "normal form initial point: " << normalFormInitialPoint << endl;
    
    vector<double> nfX, nfY;
    double nfDt = 0.1;
    CVector normalFormInitialPointC({normalFormInitialPoint[0], normalFormInitialPoint[1], normalFormInitialPoint[2], normalFormInitialPoint[3]});
    auto nfCoordsInitialPoint = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(normalFormInitialPointC), 100);

    for(double t = 0; t <= 2*timeLeft; t += nfDt)
    {
        CVector normalFormSolC = normalForm.solution(t, nfCoordsInitialPoint);
        normalFormSolC = testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormSolC));
        LDVector normalFormSol({normalFormSolC[0].real(), normalFormSolC[1].real()});

        nfX.push_back(normalFormSol[0]);
        nfY.push_back(normalFormSol[1]);
    }
    
    // trajectory from point symmetric to the one calculated by newtons method
    LDVector symmetricPoint = {-normalFormInitialPoint[0], normalFormInitialPoint[1], normalFormInitialPoint[2], -normalFormInitialPoint[3]};
    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, symmetricPoint, 10, 0, intTime);

    int N = 1000; // numer of samples in a unit of time
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

    cout << "diffLeft:\n " <<  (nfX.front() - solverX2.front()) << " " << (nfY.front() - solverY2.front()) << endl;
    cout << "diffRight:\n" << (nfX.back() - solverX.front()) << " " << (nfY.back() - solverY.front()) << endl;

    Plot2D plot, plotCloseUp;
    initializePlots(plot, plotCloseUp, normalFormInitialPoint);

    addCurveToPlots(plot, plotCloseUp, solverX, solverY, "integrated solution");
    addCurveToPlots(plot, plotCloseUp, solverX2, solverY2, "symmetric trajectory");
    addCurveToPlots(plot, plotCloseUp, nfX, nfY, "normal form solution");

    string title = "period: " + to_string(period) + ", int time: " + to_string(intTime) + ", epsilon: " + to_string(eps);
    showAndSavePlots(plot, plotCloseUp, title);

    return 0;
}