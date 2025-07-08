#pragma once

#include "capd/capdlib.h"
#include "../../../normal-forms/source/typedefs.h"
#include "../../../normal-forms/source/containers/Polynomial.hpp"
#include "../../../normal-forms/source/NormalFormFinder/PseudoNormalForm.hpp"
#include "../../../normal-forms/source/Diagonalization/Diagonalization.hpp"


// Newton's method inverse calculation
CVector inverseNewton(CVector initialPoint, const Polynomial<capd::Complex> &phi, const CVector &x);
// phi = id - Q
// phi^{-1} = id + Q + Q^2 + Q^3 + ... + Q^deg
CVector inverse(const Polynomial<capd::Complex> &phi, const CVector &x);

// calculating derivative of inverse of phi in x as an inverse of derivative in x
CMatrix inverseDer(const Polynomial<capd::Complex> &phi, const CVector &x);

// transformations from/to the space where system is in the normal form
CVector toNF(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diag);
CVector fromNF(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diag);

// derivatives
CMatrix toNF_der(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diag);
CMatrix fromNF_der(const CVector &X, const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diag);

// derivative of fromNF(normalForm.solution(t, toNF(point)))
capd::LDMatrix normalFormOriginalSolutionDerivative(double t, const capd::LDVector &point, const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diag);

// derivative of normalForm.solution(t, point)
CMatrix normalFormSolutionDerivative(double t, CVector point, const PseudoNormalForm &normalForm);