#include "test_cases_collection.h"

#include "../../../normal-forms/source/typedefs.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

using namespace std;
using namespace capd;
using capd::autodiff::Node;

void TestCasesCollection::generate_diagonal_matrix_test()
{
    diagonal_matrix.f = CMap(
        "par: p1, p2, p3, p4, a1, a2, a3;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "p1*x1 + a1*x1*x3*x4,"
            "p2*x2 + a2*x2*x4,"
            "p3*x3 + a3*x2*x2,"
            "p4*x4 + x1*x1*x1 + x2*x3;", 
            maxDerivative);

    diagonal_matrix.f.setParameter("p1", Complex(1, 1));
    diagonal_matrix.f.setParameter("p2", Complex(-1, -1));
    diagonal_matrix.f.setParameter("p3", Complex(1, -1));
    diagonal_matrix.f.setParameter("p4", Complex(-1, 1));
    diagonal_matrix.f.setParameter("a1", 2);
    diagonal_matrix.f.setParameter("a2", 5);
    diagonal_matrix.f.setParameter("a3", 1i);
    diagonal_matrix.f.setDegree(maxDerivative);

    diagonal_matrix.p = CVector({0, 0, 0, 0});
    diagonal_matrix.lambda = CMatrix({ {Complex(1, 1), 0, 0, 0}, {0, Complex(-1, -1), 0, 0}, {0, 0, Complex(1, -1), 0}, {0, 0, 0, Complex(-1, 1)} });
    diagonal_matrix.J = CMatrix({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });
    diagonal_matrix.invJ = CMatrix({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });
}

void pcr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
    Node mu = params[0]; // relative mass of Jupiter
    Node mj = 1 - mu; // relative mass of the Sun

    // in[0] = x, in[1] = vx, in[2] = y, in[3] = vy
    Node xMu = in[0] + mu;
    Node xMj = in[0] - mj;
    Node ySquare = in[2]*in[2];
    Node xMuSquareYSquare = xMu*xMu + ySquare;
    Node xMjSquareYSquare = xMj*xMu + ySquare;
    Node factor1 = mj / ( xMuSquareYSquare * sqrt(xMuSquareYSquare) );
    Node factor2 = mu / ( xMjSquareYSquare * sqrt(xMjSquareYSquare) );

    out[0] = in[1];
    out[1] = in[0] - xMu*factor1 - xMj*factor2 + 2*in[3];
    out[2] = in[3];
    out[3] = in[2]*(1 - factor1 - factor2) - 2*in[1];
}

void TestCasesCollection::generate_PCR3BP_test()
{
    int dim=4, noParam=1;
    PCR3BP.f = CMap(pcr3bpVectorField, dim, dim, noParam, maxDerivative);
    PCR3BP.f.setParameter(0, 0.5); // mu parameter
    PCR3BP.f.setDegree(maxDerivative);

    PCR3BP.p = CVector({0, 0, 0.866025403784439, 0}); // L4

    PCR3BP.lambda = CMatrix({ 
        {Complex(0.632075, 0.94843), 0, 0, 0}, 
        {0, Complex(-0.632075, -0.94843), 0, 0}, 
        {0, 0, Complex(0.632075, -0.94843), 0}, 
        {0, 0, 0, Complex(-0.632075, 0.94843)} });

    PCR3BP.J = CMatrix({ 
        {Complex(-0.677045, -0.120902), Complex(0, 1.56774), Complex(-0.417702, -0.958065), Complex(-0.660842, 0.440414)},
        {Complex(0.677045, 0.120902), Complex(0, 1.56774), Complex(-0.417702, -0.958065), Complex(0.660842, -0.440414)},
        {Complex(-0.677045, 0.120902), Complex(0, -1.56774), Complex(-0.417702, 0.958065), Complex(-0.660842, -0.440414)},
        {Complex(0.677045, -0.120902), Complex(0, -1.56774), Complex(-0.417702, 0.958065), Complex(0.660842, 0.440414)} 
    });

    PCR3BP.invJ = CMatrix({ 
        {Complex(-0.29122, 0.436975), Complex(0.29122, -0.436975), Complex(-0.29122, -0.436975), Complex(0.29122, 0.436975)},
        {Complex(-0.365758, -0.159465), Complex(-0.365758, -0.159465), Complex(-0.365758, 0.159465), Complex(-0.365758, 0.159465)},
        {-0.598513, -0.598513, -0.598513, -0.598513},
        {Complex(-0.0799453, -0.44769), Complex(0.0799453, 0.44769), Complex(-0.0799453, 0.44769), Complex(0.0799453, -0.44769)}
    });
}
