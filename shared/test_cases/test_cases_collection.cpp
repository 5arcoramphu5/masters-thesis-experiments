#include "test_cases_collection.h"

#include "../../../normal-forms/source/typedefs.h"

using namespace std;
using namespace capd;
using capd::autodiff::Node;

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

    diagonal_matrix = TestCase("diagonal_matrix", diagonalVectorField, 5, p, J, invJ, lambda, maxDerivative);
    diagonal_matrix.diagonalization.setParameter(0, Complex(1, 1));
    diagonal_matrix.diagonalization.setParameter(1, Complex(-1, -1));
    diagonal_matrix.diagonalization.setParameter(2, Complex(1, -1));
    diagonal_matrix.diagonalization.setParameter(3, Complex(-1, 1));
    diagonal_matrix.diagonalization.setParameter(4, 1i);
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
        {Complex( 0.121676201410512 , 0.564560233479915 ) ,  Complex( 1.21812018142403 , -0.499272846617592 ) ,  Complex( -0.611382974395644 , 0.629661137669460 ) ,  Complex( 0.552652795264054 , 0.373210434173016 ) ,  },
        {Complex( -0.619542857560677 , -0.298619870363536 ) ,  Complex( 0.421842688126384 , -1.50992256334870 ) ,  Complex( 0.144503592093257 , 1.03512414769062 ) ,  Complex( -0.754974145930065 , 0.246353969085843 ) ,  },
        {Complex( -0.554464141199960 , -0.114450545144107 ) ,  Complex( 0.479012706985424 , -1.19835992210033 ) ,  Complex( -0.612015114862388 , 0.604705272381974 ) ,  Complex( -0.370572528956227 , -0.538561555438970 ) ,  },
        {Complex( -0.407756162269007 , -0.865173057898415 ) ,  Complex( -2.10485198676731 , 0.568316665220325 ) ,  Complex( 1.43771793515155 , 0.213501788479026 ) ,  Complex( 0.351740568807508 , -1.04689919105447 ) ,  }});

    CMatrix invJ({
        {Complex( -0.349978887198728 , -0.518251347185310 ) ,  Complex( -0.162899282799330 , 0.499219669000705 ) ,  Complex( -0.525526721991652 , -0.361603542688430 ) ,  Complex( -0.357940656558002 , 0.120262056951474 ) ,  },
        {Complex( 0.340905753627265 , -0.331009746627618 ) ,  Complex( 0.395176782768370 , 0.0551667785446674 ) ,  Complex( -0.340676818786865 , 0.344795013212806 ) ,  Complex( -0.0421450968509044 , -0.283804468585034 ) ,  },
        {Complex( 0.270312039062428 , -0.659504221578514 ) ,  Complex( 0.576439398254525 , -0.161045838511314 ) ,  Complex( -0.675127975012814 , 0.269864564817615 ) ,  Complex( 0.112185293942147 , -0.415496242308252 ) ,  },
        {Complex( 0.529417572978029 , 0.114102119525623 ) ,  Complex( -0.197459626456834 , -0.409666982572198 ) ,  Complex( 0.111680492623941 , 0.541044416639114 ) ,  Complex( 0.295807480822040 , 0.139414099929843 ) ,  }});

    PCR3BP_L4 = TestCase("PCR3BP_L4", pcr3bpVectorField, 1, p, J, invJ, lambda, maxDerivative);
    PCR3BP_L4.diagonalization.setParameter(0, 0.5); // mu parameter
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

    PCR3BP_L1 = TestCase("PCR3BP_L1", pcr3bpVectorField, 1, p, J, invJ, lambda, maxDerivative);
    PCR3BP_L1.diagonalization.setParameter(0, 0.5); // mu parameter
}

TestCasesCollection::TestCasesCollection(int maxDerivative): maxDerivative(maxDerivative)
{
    generate_diagonal_matrix_test();
    generate_PCR3BP_L1_test();
    generate_PCR3BP_L4_test();

    PCR3BP_dmap = LDMap(pcr3bpVectorField, 4, 4, 1, maxDerivative);
    PCR3BP_dmap.setParameter(0, 0.5); // mu parameter
    PCR3BP_dmap.setDegree(maxDerivative);
}