#include <iostream>
#include <iomanip>
#include <vector>
#include <sciplot/sciplot.hpp>

#include "../../shared/calculations/calculations.h"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

// INTEGRATION --------------------------------------------------------

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &initialPoint, int order, double initTime, double finalTime)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(initTime);
    LDVector point(initialPoint);
    timeMap(finalTime, point, solution);
    return solution;
}

// PLOTTING ----------------------------------------------------------

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

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(1000, 1000);

    canvas.show();
    canvas.save("experiments/periodic-orbit-for-given-period-normal-form/" + filename);
}

void showAndSavePlots(Plot2D &plot, Plot2D &plotCloseUp)
{
    showAndSavePlot(plot, "plot.pdf");
    showAndSavePlot(plotCloseUp, "plotCloseUp.pdf");
}

// otrzymuje (y, v_x) 
// oblicza punkt końcowy dla casu halfPeriod przy punkcie początkowym (0, y, v_x, 0)
// zwraca rzutowanie punktu końcowego na zmienne (x, v_y)
// celem jest znalezienie miejsca zerowego, więc punktu (0, y', v_x', 0)
// zwraca wynik i pochodną
pair<LDVector, LDMatrix> computeF(const LDVector &initialPoint, double halfPeriod, double intTime, LDMap &map, int order, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    // compute value of the function
    LDVector point({0, initialPoint[0], initialPoint[1], 0});
    LDMatrix intDer;
    point = timeMap(intTime, point, intDer);

    CVector lastPoint({point[0], point[1], point[2], point[3]});
    auto newInitialPoint = inverse(normalForm.getPhi(), diag.toDiag(lastPoint));

    CVector normalFormSol = normalForm.solution(halfPeriod - intTime, newInitialPoint);
    normalFormSol = diag.toOriginal(normalForm.getPhi()(normalFormSol));
    LDVector result({normalFormSol[0].real(), normalFormSol[3].real()});

    // compute derivative
    auto der = normalFormOriginalSolutionDerivative(halfPeriod - intTime, point, normalForm, diag) * intDer;
    
    return make_pair(result, der);
}

LDVector getInitialPoint(double halfPeriod, double intTime, const LDVector &approximation, LDMap &map, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag, double eps)
{
    int solverOrder = 25;
    int maxIters = 20;
    LDVector point(approximation);

    auto bestPoint = point;
    optional<LDVector> bestF = std::nullopt;

    long double damp = 1;

    for(int i = 0; i < maxIters; ++i)
    {
        auto [f, f_der] = computeF(point, halfPeriod, intTime, map, solverOrder, normalForm, diag);

        cout << f << endl;

        if(!bestF.has_value() || abs(f) < abs(bestF.value()))
        {
            bestF = f;
            bestPoint = point;
        }

        if(abs(f[0]) < eps && abs(f[1]) < eps)
            break;

        // trunceting matrix to take projections into account
        LDMatrix der({ {f_der[0][1], f_der[0][2]}, {f_der[3][1], f_der[3][2]}});

        // inverse of 2x2 matrix
        long double det = der[0][0]*der[1][1] - der[0][1]*der[1][0];
        LDMatrix invDer({   {der[1][1]/det, -der[0][1]/det}, 
                            {-der[1][0]/det, der[0][0]/det} });

        // invDer * f
        LDVector invDerF = invDer * f;

        // damped Newton's method
        while(damp > 1e-10)
        {
            cout << ".";
            LDVector newPoint = point - damp*invDerF;
            auto [newF, newFDer] = computeF(newPoint, halfPeriod, intTime, map, solverOrder, normalForm, diag);

            if(abs(newF) > abs(f))
                damp *= 0.5; // if the function value is worse, reduce the step
            else
                break;
        }
        
        point = point - damp*invDerF;
    }

    return LDVector({0, bestPoint[0], bestPoint[1], 0});
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

    LDVector approximation({-0.576003306698015024546, -0.451540030090179645464}); // w przybliżeniu orbita okresowa zaczyna się w (0, x, y, 0)
    auto initialPoint = getInitialPoint(period/2, intTime, approximation, testCases.PCR3BP_dmap, normalForm, testCase.diagonalization, eps);
    cout << "initial point: " << initialPoint << endl;
    CVector initialPointC({initialPoint[0], initialPoint[1], initialPoint[2], initialPoint[3]});

    auto intSolution = integrateSolution(testCases.PCR3BP_dmap, initialPoint, 10, 0, intTime);

    int N = 1000; // numer of samples in a unit of time
    double intDt = 1. / N;

    Plot2D plot, plotCloseUp;
    initializePlots(plot, plotCloseUp, intSolution(intTime));

    vector<double> solverX, solverY, solverX2, solverY2;

    for(double t = 0; t <= intTime; t += intDt)
    {
        solverX.push_back(intSolution(t)[0]);
        solverY.push_back(intSolution(t)[1]);
        // symmetry with respect to mass
        solverX2.push_back(-intSolution(t)[0]); 
        solverY2.push_back(intSolution(t)[1]); 
    }
    addCurveToPlots(plot, plotCloseUp, solverX, solverY, "integrated solution");
    addCurveToPlots(plot, plotCloseUp, solverX2, solverY2, "symmetric trajectory");

    if(timeLeft > 0)
    {
        vector<double> nfX, nfY;
        double nfDt = timeLeft / N;
        cout << nfDt << endl;
        
        CVector lastPoint = CVector({intSolution(intTime)[0], intSolution(intTime)[1], intSolution(intTime)[2], intSolution(intTime)[3]});
        auto newInitialPoint = inverse(normalForm.getPhi(), testCase.diagonalization.toDiag(lastPoint));

        for(double t = 0; t <= 2*timeLeft; t += nfDt)
        {
            CVector normalFormSolC = normalForm.solution(t, newInitialPoint);
            normalFormSolC = testCase.diagonalization.toOriginal(normalForm.getPhi()(normalFormSolC));
            LDVector normalFormSol({normalFormSolC[0].real(), normalFormSolC[1].real()});
 
            nfX.push_back(normalFormSol[0]);
            nfY.push_back(normalFormSol[1]);
        }

        cout << "diffLeft:\n " <<  (nfX.front() - solverX.back()) << " " << (nfY.front() - solverY.back()) << endl;
        cout << "diffRight:\n" << (nfX.back() - solverX2.back()) << " " << (nfY.back() - solverY2.back()) << endl;

        addCurveToPlots(plot, plotCloseUp, nfX, nfY, "normal form solution");
    }

    showAndSavePlots(plot, plotCloseUp);

    return 0;
}