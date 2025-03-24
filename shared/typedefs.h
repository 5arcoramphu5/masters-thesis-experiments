#pragma once

#include "../../normal-forms/source/typedefs.h"
#include "capd/dynsys/lib.h"
#include "capd/poincare/lib.h"

typedef capd::dynsys::OdeSolver<CMap> COdeSolver;
typedef capd::poincare::TimeMap<COdeSolver> CTimeMap;