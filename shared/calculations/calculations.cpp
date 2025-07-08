#include "calculations.h"
#include "capd/matrixAlgorithms/lib.h"

using namespace std;
using namespace capd;
using namespace capd::matrixAlgorithms;

#define INVERSE_DEG 500
#define INVERSE_PRECISION 1e-13

// Newtons method inverse calculation
// phi(inv) = x
// F(inv) = phi(inv) - x  
// DF(inv) = phi.derivative(inv)
CVector inverseNewton(CVector initialPoint, const Polynomial<Complex> &phi, const CVector &x)
{
    CVector inv(initialPoint);

    while(true)
    {
        CVector F = phi(inv) - x;
        CMatrix DF = phi.derivative(inv);

        if(abs(F[0]) < INVERSE_PRECISION && abs(F[1]) < INVERSE_PRECISION && abs(F[2]) < INVERSE_PRECISION && abs(F[3]) < INVERSE_PRECISION)
            break;

        // invDerF = (f_der)^{-1} * f
        // f_der * invDerF = f - rownanie liniowe
        CVector invderF = matrixAlgorithms::gauss(DF, F);

        inv -= invderF;
    }

    return inv;
}

// phi = id - Q
// Q = id - phi
// phi^{-1} = id + Q + Q^2 + Q^3 + ... + Q^deg
CVector inverse(const Polynomial<Complex> &phi, const CVector &x)
{
    CVector result = x;
    CVector curr = x;
    for(int i = 0; i < INVERSE_DEG; ++i)
    {
        curr = curr - phi(curr); // Q(curr)
        result += curr;
    }

    return inverseNewton(result, phi, x); // correct the result with Newton's method
}

// calculating derivative of inverse of phi in x as an inverse of derivative in x
CMatrix inverseDer(const Polynomial<Complex> &phi, const CVector &x)
{
    CMatrix der = phi.derivative(inverse(phi, x));
    CMatrix derInv = gaussInverseMatrix(der);
    return derInv;
}

// derivative of normalForm.solution(t, point)
CMatrix normalFormSolutionDerivative(double t, CVector point, const PseudoNormalForm &normalForm)
{
    CVector mulV({ point[0]*point[1], point[2]*point[3]});

    auto a1 = normalForm.geta1()(mulV);
    auto a2 = normalForm.geta2()(mulV); 

    auto exp1 = exp(a1[0]*t);
    auto exp2 = exp(-a1[0]*t);
    auto exp3 = exp(a2[0]*t);
    auto exp4 = exp(-a2[0]*t);

    auto derA1 = normalForm.geta1().derivative(mulV);
    auto derA2 = normalForm.geta2().derivative(mulV);

    CVector derA1mul({
        derA1[0][0]*point[1],
        derA1[0][0]*point[0],
        derA1[0][1]*point[3],
        derA1[0][1]*point[2]
    });

    CVector derA2mul({
        derA2[0][0]*point[1],
        derA2[0][0]*point[0],
        derA2[0][1]*point[3],
        derA2[0][1]*point[2]
    });

    Complex one = 1;
    CMatrix result({
        { (one + derA1mul[0]*point[0]*t) * exp1, point[0]*exp1*derA1mul[1]*t, point[0]*exp1*derA1mul[2]*t, point[0]*exp1*derA1mul[3]*t },
        { -point[1]*exp2*derA1mul[0]*t, (one - derA1mul[1]*point[1]*t) * exp2, -point[1]*exp2*derA1mul[2]*t, -point[1]*exp2*derA1mul[3]*t },
        { point[2]*exp3*derA2mul[0]*t, point[2]*exp3*derA2mul[1]*t, (one + derA2mul[2]*point[2]*t) * exp3, point[2]*exp3*derA2mul[3]*t },
        { -point[3]*exp4*derA2mul[0]*t, -point[3]*exp4*derA2mul[1]*t, -point[3]*exp4*derA2mul[2]*t, (one - derA2mul[3]*point[3]*t) * exp4}
    });

    // cout << "normalFormSolutionDerivative: \n" << result << endl;
    return result;
}

// phi^{-1}(diag.toDiag(X))
CVector toNF(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    CVector diagPoint = diag.toDiag(X);
    CVector invPhiDiag = inverse(normalForm.getPhi(), diagPoint);
    return normalForm.solution(0, invPhiDiag);
}

// diag.toOrig(phi(X))
CVector fromNF(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    CVector phiX = normalForm.getPhi()(X);
    return diag.toOriginal(phiX);
}

// derivative of phi^{-1}(diag.toDiag(X))
CMatrix toNF_der(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    return inverseDer(normalForm.getPhi(), diag.toDiag(X)) * diag.getJ();
}

// derivative of diag.toOrig(phi(X))
CMatrix fromNF_der(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    return diag.getinvJ() * normalForm.getPhi().derivative(X);
}

// derivative of fromNF(normalForm.solution(t, toNF(point)))
LDMatrix normalFormOriginalSolutionDerivative(double t, const LDVector &point, const PseudoNormalForm &normalForm, const Diagonalization<Complex> &diag)
{
    CVector pointC({point[0], point[1], point[2], point[3]});

    CVector nfCoords = toNF(pointC, normalForm, diag);
    CVector nfSol = normalForm.solution(t, nfCoords);

    CMatrix resultComplex = fromNF_der(nfSol, normalForm, diag) * normalFormSolutionDerivative(t, nfCoords, normalForm) * toNF_der(pointC, normalForm, diag);

    LDMatrix result(4, 4);
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            result[i][j] = resultComplex[i][j].real();

    return result;
}