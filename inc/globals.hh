#ifndef _H_GLOBALS_
#define _H_GLOBALS_
#include "types.hh"

namespace RS {

class Globals {
private:
  Globals() {}

public:
  static real EPS;
  static real ZERO;
  static real SMALL;

  static real SEARCHPREC; // threshold for recursive search in spatial distance

  static real SPACEEQUAL;
  static real TIMEEQUAL;
  static real T0EQUAL;
  static real TAUEQUAL;
  static real TAUMIN;
  static real DETMIN;
  static real RECPOINTEQUAL;   // minimum distance in R4
  static real RCPCLUSTERMERGE; // merge distance in R4
};
} // namespace FT

#endif //_H_GLOBALS_
