#ifndef _H_JOBPARAMS_
#define _H_JOBPARAMS_

#include "globals.hh"
#include "types.hh"
#include <sstream>

namespace RS {
//-----------------------------------------------------------------------------------------------//
struct DataParams {
  DataParams(){};
  DataParams(const Vec3r&              dMin,
             const Vec3r&              dMax,
             const Vec3i               samples,
             const char*               flow,
             const char*               savefile,
             std::vector<std::string>* p_loadFiles = nullptr,
             Vec2r*                    p_timeRange = nullptr)
  : dMin(dMin)
  , dMax(dMax)
  , samples(samples)
  , flowName(flow)
  , saveFilename(savefile) {
    if (nullptr != p_loadFiles)
      loadFiles = *p_loadFiles;
    if (nullptr != p_timeRange)
      timeRange = *p_timeRange;
  };

  Vec3r       dMin         = Vec3r(-1.0, -1.0, -1.0); // domainMin
  Vec3r       dMax         = Vec3r(1.0, 1.0, 1.0);    // domainMax
  Vec3i       samples      = Vec3i(0, 0, 0);
  std::string flowName     = "none";
  std::string saveFilename = "savefile.ply";

  std::vector<std::string> loadFiles;
  Vec2r                    timeRange = Vec2r(0.0, 0.0);

  virtual std::string toStr(void) const {
    std::stringstream os;
    os << "domain: [" << dMin[0] << ", " << dMax[0] << "]x";
    os << "[" << dMin[1] << ", " << dMax[1] << "]x";
    os << "[" << dMin[2] << ", " << dMax[2] << "]\n";
    os << "samples=[" << samples[0] << "x" << samples[1] << "x" << samples[2]
       << "]\n";
    os << "flow: " << flowName << "\n";
    os << "saveFile: " << saveFilename << "\n";
    return os.str();
  }
};

//-----------------------------------------------------------------------------------------------//
struct SearchParams {
  SearchParams(const real& t0_min,
               const real& t0_max,
               const real& tau_min,
               const real& tau_max,
               const real& dt,
               const real& lowPrec,
               const real& highPrec)
  : t0_min(t0_min)
  , t0_max(t0_max)
  , tau_min(tau_min)
  , tau_max(tau_max)
  , dt(dt)
  , lowPrec(lowPrec)
  , highPrec(highPrec){};

  real t0_min   = 0.0;
  real t0_max   = 14.8;
  real tau_min  = 0.2;
  real tau_max  = 15.0;
  real dt       = 0.2;
  real lowPrec  = Globals::SEARCHPREC * 100.0;
  real highPrec = Globals::SEARCHPREC;

  virtual std::string toStr(void) const {
    std::stringstream os;
    os << "t0_min=" << t0_min << " t0_max=" << t0_max << " tau_min=" << tau_min
       << " tau_max=" << tau_max << " dt=" << dt;
    return os.str();
  }
};
}

#endif //_H_JOBPARAMS_
