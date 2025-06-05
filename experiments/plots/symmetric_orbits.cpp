#include <iostream>
#include <sciplot/sciplot.hpp>
#include "capd/capdlib.h"
#include "../../shared/test_cases/test_cases_collection.h"

using namespace std;
using namespace capd;
using namespace sciplot;

LDTimeMap::SolutionCurve integrateSolution(LDMap &map, LDVector &point, int order, double finalTime)
{
    LDOdeSolver solver(map, order);
    LDTimeMap timeMap(solver);

    LDTimeMap::SolutionCurve solution(0);
    timeMap(finalTime, point, solution);
    return solution;
}


Plot2D initializePlot()
{
    Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.palette("set1");

    plot.xrange(-1.25, 1.25);
    plot.yrange(-0.9, 0.9);
    plot.size(1000, 1000);

    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // drawing libration points
    vector<double> librationPointsX({0, 0, 0, 1.19841, -1.19841});
    vector<double> librationPointsY({0, 0.866025403784438646763723170753, -0.866025403784438646763723170753, 0, 0});
    plot.drawPoints(librationPointsX, librationPointsY).label("punkty libracji");
    vector<double> primatesX({-0.5, 0.5});
    vector<double> primatesY({0, 0});
    plot.drawPoints(primatesX, primatesY).label("Słońce i Jowisz");

    return plot;
}

void showAndSavePlot(Plot2D &plot, string filename)
{
    Figure fig({{plot}});
    Canvas canvas({{fig}});
    canvas.size(500, 500);

    canvas.show();
    canvas.save("experiments/plots/symmetric_orbits/" + filename);
}

void drawOrbitAndTimeSymmetric(LDVector &point, Plot2D &plot, TestCasesCollection &testCases)
{
    double finalTime = 3;
    auto solution = integrateSolution(testCases.PCR3BP_dmap, point, 10, finalTime);

    auto lastPoint = solution(finalTime);
    LDVector symmetricPoint({lastPoint[0], -lastPoint[1], -lastPoint[2], lastPoint[3]});
    auto symmetricSolution = integrateSolution(testCases.PCR3BP_dmap, symmetricPoint, 10, finalTime);

    int N = 1000; // number of samples in a unit of time
    double dt = 1./N;
    vector<double> x, y, symX, symY;
    x.reserve(N*finalTime);
    y.reserve(N*finalTime);
    symX.reserve(N*finalTime);
    symY.reserve(N*finalTime);

    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(solution(t)[0]);
        y.push_back(solution(t)[1]);
        symX.push_back(symmetricSolution(t)[0]);
        symY.push_back(symmetricSolution(t)[1]);
    }
    
    plot.drawCurve(x, y).label("orbita");
    plot.drawCurve(symX, symY).label("orbita symetryczna");
}

void drawTimeSymmetricOrbit(LDVector &point, Plot2D &plot, TestCasesCollection &testCases)
{
    double finalTime = 3;
    auto solution = integrateSolution(testCases.PCR3BP_dmap, point, 10, finalTime);

    auto lastPoint = solution(finalTime);
    LDVector symmetricPoint({lastPoint[0], -lastPoint[1], -lastPoint[2], lastPoint[3]});
    auto symmetricSolution = integrateSolution(testCases.PCR3BP_dmap, symmetricPoint, 10, finalTime);

    int N = 1000; // number of samples in a unit of time
    double dt = 1./N;
    vector<double> x, y, symX, symY;
    x.reserve(N*finalTime*2);
    y.reserve(N*finalTime*2);

    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(symmetricSolution(t)[0]);
        y.push_back(symmetricSolution(t)[1]);
    }
    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(solution(t)[0]);
        y.push_back(solution(t)[1]);
    }

    plot.drawCurve(x, y).label("orbita symetryczna");
}

void drawOrbitAndMassSymmetric(LDVector &point, Plot2D &plot, TestCasesCollection &testCases)
{
    double finalTime = 5;
    auto solution = integrateSolution(testCases.PCR3BP_dmap, point, 10, finalTime);

    auto lastPoint = solution(finalTime);
    LDVector symmetricPoint({-lastPoint[0], lastPoint[1], lastPoint[2], -lastPoint[3]});
    auto symmetricSolution = integrateSolution(testCases.PCR3BP_dmap, symmetricPoint, 10, finalTime);

    int N = 1000; // number of samples in a unit of time
    double dt = 1./N;
    vector<double> x, y, symX, symY;
    x.reserve(N*finalTime);
    y.reserve(N*finalTime);
    symX.reserve(N*finalTime);
    symY.reserve(N*finalTime);

    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(solution(t)[0]);
        y.push_back(solution(t)[1]);
        symX.push_back(symmetricSolution(t)[0]);
        symY.push_back(symmetricSolution(t)[1]);
    }
    
    plot.drawCurve(x, y).label("orbita");
    plot.drawCurve(symX, symY).label("orbita symetryczna");
}

void drawMassSymmetricOrbit(LDVector &point, Plot2D &plot, TestCasesCollection &testCases)
{
    double finalTime = 5;
    auto solution = integrateSolution(testCases.PCR3BP_dmap, point, 10, finalTime);

    auto lastPoint = solution(finalTime);
    LDVector symmetricPoint({-lastPoint[0], lastPoint[1], lastPoint[2], -lastPoint[3]});
    auto symmetricSolution = integrateSolution(testCases.PCR3BP_dmap, symmetricPoint, 10, finalTime);

    int N = 1000; // number of samples in a unit of time
    double dt = 1./N;
    vector<double> x, y, symX, symY;
    x.reserve(N*finalTime*2);
    y.reserve(N*finalTime*2);

    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(symmetricSolution(t)[0]);
        y.push_back(symmetricSolution(t)[1]);
    }
    for(double t = 0; t <= finalTime; t += dt)
    {
        x.push_back(solution(t)[0]);
        y.push_back(solution(t)[1]);
    }

    plot.drawCurve(x, y).label("orbita symetryczna");
}

int main()
{
    TestCasesCollection testCases(10);

    // para orbit wykazujących symetrię względem odwrocenia czasu
    Plot2D plot1 = initializePlot();
    LDVector point1({0.5, -0.57, -0.3, -0.1});
    drawOrbitAndTimeSymmetric(point1, plot1, testCases);
    showAndSavePlot(plot1, "symmetric_orbit_time_1.pdf");

    // orbita symetryczna względem odwrocenia czasu
    Plot2D plot2 = initializePlot();
    LDVector point2({0.1, 0, 0, 0.1});
    drawTimeSymmetricOrbit(point2, plot2, testCases);
    plot2.drawPoints(vector<double>({0.1}), vector<double>({0})).label("punkt (0.1, 0, 0, 0.1)").pointType(7);
    showAndSavePlot(plot2, "symmetric_orbit_time_2.pdf");

    // para orbit symetrycznych względem masy
    Plot2D plot3 = initializePlot();
    LDVector point3({0.5, -0.57, -0.3, -0.1});
    drawOrbitAndMassSymmetric(point3, plot3, testCases);
    showAndSavePlot(plot3, "symmetric_orbit_mass_1.pdf");

    // orbita symetryczna względem masy
    Plot2D plot4 = initializePlot();
    LDVector point4({0, 0.2, 0.2, 0});
    drawMassSymmetricOrbit(point4, plot4, testCases);
    plot4.drawPoints(vector<double>({0}), vector<double>({0.2})).label("punkt (0, 0.2, 0.2, 0)").pointType(7);
    showAndSavePlot(plot4, "symmetric_orbit_mass_2.pdf");

    return 0;
}