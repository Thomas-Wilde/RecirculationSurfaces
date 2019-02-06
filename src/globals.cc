#include "globals.hh"

using namespace RS;

real Globals::EPS             = std::numeric_limits<double>::epsilon();
real Globals::ZERO            = 1.0e3 * Globals::EPS;
real Globals::SMALL           = 1.0e6 * Globals::EPS;
real Globals::SEARCHPREC      = 1.0e-6;
real Globals::SPACEEQUAL      = 5.0e-5;
real Globals::TIMEEQUAL       = 5.0e-5;
real Globals::T0EQUAL         = 5.0e-5;
real Globals::TAUEQUAL        = 5.0e-5;
real Globals::TAUMIN          = 1.0e-3;
real Globals::DETMIN          = 1.0e-6;
real Globals::RECPOINTEQUAL   = 5.0e-5;
real Globals::RCPCLUSTERMERGE = 5.0e-3;