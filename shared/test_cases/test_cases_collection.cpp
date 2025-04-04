#include "test_cases_collection.h"

#include "../../../normal-forms/source/typedefs.h"
#include "../../../normal-forms/source/NormalFormFinder/NormalFormFinder.hpp"

using namespace std;
using namespace capd;
using capd::autodiff::Node;

TestCasesCollection::TestCasesCollection(int maxDerivative): maxDerivative(maxDerivative)
{
    generate_diagonal_matrix_test();
    generate_PCR3BP_L1_test();
    generate_PCR3BP_L4_test();
}

void diagonalVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
    out[0] = params[0]*in[0] + 2*in[0]*in[2]*in[3];
    out[1] = params[1]*in[1] + 5*in[1]*in[3];
    out[2] = params[2]*in[2] + params[4]*in[1]*in[1];
    out[3] = params[3]*in[3] + in[0]*in[0]*in[0] + in[1]*in[2];
}

void TestCasesCollection::generate_diagonal_matrix_test()
{
    CVector p({0, 0, 0, 0});
    CMatrix lambda({ {Complex(1, 1), 0, 0, 0}, {0, Complex(-1, -1), 0, 0}, {0, 0, Complex(1, -1), 0}, {0, 0, 0, Complex(-1, 1)} });
    CMatrix J = CMatrix::Identity(4);
    CMatrix invJ = CMatrix::Identity(4);

    Diagonalization<Complex> diagonalization(diagonalVectorField, 5, p, J, invJ, lambda, maxDerivative);
    diagonalization.setParameter(0, Complex(1, 1));
    diagonalization.setParameter(1, Complex(-1, -1));
    diagonalization.setParameter(2, Complex(1, -1));
    diagonalization.setParameter(3, Complex(-1, 1));
    diagonalization.setParameter(4, 1i);

    diagonal_matrix = TestCase("diagonal_matrix", diagonalization);
}

void pcr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{       
    Node mu = params[0]; // relative mass of Jupiter
    Node mj = 1 - mu; // relative mass of the Sun

    Node xMu = in[0] + mu;
    Node xMj = in[0] - mj;
    Node xMuSquare = xMu^2;
    Node xMjSquare = xMj^2;
    Node ySquare = in[1]^2;
    Node factor1 = mj / ((xMuSquare + ySquare)^1.5);
    Node factor2 = mu / ((xMjSquare + ySquare)^1.5);

    out[0] = in[2];
    out[1] = in[3];
    out[2] = in[0] - xMu*factor1 - xMj*factor2 + 2*in[3];
    out[3] = in[1] * (1 - factor1 - factor2) - 2*in[2];
}

void TestCasesCollection::generate_PCR3BP_L4_test()
{
    CVector p({0, 0.866025403784438646763723170753, 0,  0}); // L4

    CMatrix lambda({ 
        {Complex(0.632075, 0.94843), 0, 0, 0}, 
        {0, Complex(-0.632075, -0.94843), 0, 0}, 
        {0, 0, Complex(0.632075, -0.94843), 0}, 
        {0, 0, 0, Complex(-0.632075, 0.94843)} });

    CMatrix J({
        {Complex(-0.677045, -0.120902), Complex(0.0000000000000000734646, 1.56774), Complex(-0.417702, -0.958065), Complex(-0.660842, 0.440414)}, 
        {Complex(0.677045, 0.120902), Complex(-0.000000000000000118226, 1.56774), Complex(-0.417702, -0.958065), Complex(0.660842, -0.440414)}, 
        {Complex(-0.677045, 0.120902), Complex(0.000000000000000203792, -1.56774), Complex(-0.417702, 0.958065), Complex(-0.660842, -0.440414)}, 
        {Complex(0.677045, -0.120902), Complex(0.0000000000000000742892, -1.56774), Complex(-0.417702, 0.958065), Complex(0.660842, 0.440414)}});

    CMatrix invJ({
        {Complex(-0.29122, 0.436975), Complex(0.29122, -0.436975), Complex(-0.29122, -0.436975), Complex(0.29122, 0.436975)}, 
        {Complex(-0.365758, -0.159465), Complex(-0.365758, -0.159465), Complex(-0.365758, 0.159465), Complex(-0.365758, 0.159465)}, 
        {Complex(-0.598513, 0.), Complex(-0.598513, 0.), Complex(-0.598513, 0.), Complex(-0.598513, 0.)}, 
        {Complex(-0.0799453, -0.44769), Complex(0.0799453, 0.44769), Complex(-0.0799453, 0.44769), Complex(0.0799453, -0.44769)}});

    Diagonalization<Complex> diagonalization(pcr3bpVectorField, 1, p, J, invJ, lambda, maxDerivative);
    diagonalization.setParameter(0, 0.5); // mu parameter

    PCR3BP_L4 = TestCase("PCR3BP_L4", diagonalization);

    PCR3BP_dmap = DMap(pcr3bpVectorField, 4, 4, 1, maxDerivative);
    PCR3BP_dmap.setDegree(maxDerivative);
}

void TestCasesCollection::generate_PCR3BP_L1_test()
{
    CVector p({0, 0, 0,  0}); // L1

    CMatrix lambda({ 
        {Complex(-3.78335, 0), 0, 0, 0}, 
        {0, Complex(3.78335, 0), 0, 0}, 
        {0, 0, Complex(0, 2.88335), 0}, 
        {0, 0, 0, Complex(0, -2.88335)} });

    CMatrix J({
        {Complex(-2.32277, 0.0000000000000000650065), Complex(-0.33955, -0.0000000000000000154575), Complex(0.516933, 0.000000000000000015734), Complex(-0.183519, -0.00000000000000000000623407)}, 
        {Complex(-0.815576, 0.0000000000000000586703), Complex(-0.000000000000000000411776, 1.47415), Complex(0.00000000000000012026, 0.138329), Complex(-0.607213, 0.0000000000000000112401)}, 
        {Complex(-2.32277, -0.00000000000000000834414), Complex(0.33955, -0.00000000000000000445094), Complex(-0.516933, -0.0000000000000000017932), Complex(-0.183519, -0.00000000000000000233778)}, 
        {Complex(-0.815576, -0.0000000000000000499585), Complex(0.0000000000000000348333, -1.47415), Complex(0.00000000000000000618679, -0.138329), Complex(-0.607213, -0.0000000000000000307705)}});

    CMatrix invJ({
        {Complex(-0.240815, 0.), Complex(0.0727819, -0.000000000000000022347), Complex(-0.240815, 0.), Complex(0.0727819, 0.000000000000000022347)},
        {Complex(-0.085493, 0.), Complex(-0.00000000000000000199701, -0.319486), Complex(0.085493, 0.), Complex(-0.00000000000000000199701, 0.319486)}, 
        {Complex(0.911087, 0.), Complex(0.00000000000000014884, -0.209856), Complex(-0.911087, 0.), Complex(0.00000000000000014884, 0.209856)}, 
        {Complex(0.32345, 0.), Complex(-0.92119, 0.), Complex(0.32345, 0.), Complex(-0.92119, 0.)}});

    Diagonalization<Complex> diagonalization(pcr3bpVectorField, 1, p, J, invJ, lambda, maxDerivative);
    diagonalization.setParameter(0, 0.5); // mu parameter

    PCR3BP_L1 = TestCase("PCR3BP_L1", diagonalization);
}