#include "jobworker.hh"
// #include "omp.h"
// #include <chrono>
// #include <thread>

using namespace std;
namespace RS {
//--------------------------------------------------------------------------//
JobWorker::JobWorker(void) {}

//--------------------------------------------------------------------------//
JobWorker::~JobWorker(void) {}

//--------------------------------------------------------------------------//
/**Uses a grid of axis parallel lines to search Recirculation Points.
Computes the recirculation points and saves them into the file, fiven in 'data'.
\b Attention: for recursive precision \c search.highPrec is used.
\param [in] data - contains information about the dataset, domain, sampling and
output files \param [in] search - contains information about the considered time
intervals and search parameters \param [in] mainDir - axis parallel direction,
in which the grid will be 'moved'; has to be one of the following: \li (1,0,0) -
for grid in yz-plane, search in x-direction \li (0,1,0) - for grid in xz-plane,
search in y-direction \li (0,0,1) - for grid in xy-plane, search in z-direction
\param [out] p_recPoints - optional; found points are stored here
*/
void
JobWorker::searchDomainByBaseGrid(const DataParams&   data,
                                  const SearchParams& search,
                                  const Vec3i&        mainDir,
                                  Flow3D*             p_flow,
                                  vector<RecPoint>*   p_recPoints) const {
  // generate the base grid, search in main-direction
  vector<Vec3r> baseGrid;
  generateBaseGrid(&baseGrid, data, mainDir);

  vector<RecPoint> recPoints =
    searchByBaseGrid(baseGrid, data, search, mainDir, search.highPrec, p_flow);

  // save result
  string             l0Filename = data.saveFilename + std::string("-hits.ply");
  RS::IORecPointFile rcpFile;
  if (0 != recPoints.size())
    rcpFile.saveFile(l0Filename, recPoints, std::string("Thomas Wilde"));

  if (p_recPoints)
    *p_recPoints = recPoints;
}

//-----------------------------------------------------------------------------------------------//
/**Computes a vector based on the sampling and the domain, which contains the
distance between two samples in x-, y- and z-direction. \param [in] data -
struct containing the sampling and domain info \return The vector containinig
the distance in x-, y- and z-direction*/
Vec3r
JobWorker::computeOffsetVector(const DataParams& data) const {
  // compute the space between samples based on domain and sample count
  Vec3r dim = data.dMax - data.dMin; // dimension of the domain
  Vec3r off(0.0, 0.0, 0.0); // offset between two samples in x,y,z-direction
  Vec3i sam = data.samples; // samples in x,y,z-direction
  for (int i = 0; i < 3; i++)
    off[i] = sam[i] > 1 ? (dim[i]) / static_cast<real>(sam[i] - 1) : 0.0;

  return off;
}

//-----------------------------------------------------------------------------------------------//
/**
\param[out] p_baseGrid - vector in which the positions of the base grid are
stored \param[in]  data - holds information about the domain and the sampling
\param[in]  searchDir - axis parallel direction, in which the grid will be
'moved'; has to be one of the following: \li (1,0,0) - for grid in yz-plane ?
search in x-direction \li (0,1,0) - for grid in xz-plane ? search in y-direction
                        \li (0,0,1) - for grid in xy-plane ? search in
z-direction \return false - if the searchDir is defined wrong
*/
bool
JobWorker::generateBaseGrid(vector<Vec3r>*    p_baseGrid,
                            const DataParams& data,
                            const Vec3i&      searchDir) const {
  // check if the search direction is given in a proper format
  if (!(searchDir[0] == 1 && searchDir[1] == 0 && searchDir[2] == 0) &&
      !(searchDir[0] == 0 && searchDir[1] == 1 && searchDir[2] == 0) &&
      !(searchDir[0] == 0 && searchDir[1] == 0 && searchDir[2] == 1)) {
    assert(false);
    return false;
  }

  Vec3i sam = data.samples; // samples in x,y,z-direction
  Vec3r off =
    computeOffsetVector(data); // offset between two samples in x,y,z-direction

  Vec3i helpVec =
    Vec3i(1, 1, 1) - searchDir;    // vector which holds the 0-dimension
  sam = sam * helpVec + searchDir; // set the samples in search direction to 1

  // create a vector with all positions on the baseGrid
  std::vector<Vec3r>& baseGrid = *p_baseGrid;
  int                 sx       = sam[0];
  int                 sy       = sam[1];
  int                 sz       = sam[2];
  Vec3r               ori      = data.dMin; // position of the first sample
  baseGrid.resize(sx * sy * sz);
  int index = 0;
  for (int i = 0; i < sx; i++)
    for (int j = 0; j < sy; j++)
      for (int k = 0; k < sz; k++)
        baseGrid[index++] = ori + (off * Vec3r(i, j, k));

  return true;
}

//--------------------------------------------------------------------------//
/**Search along consecutive lines. The lines have their origin at the points
provided by \c `baseGrid`. The search direction is provided by mainDir. The
distance between two HyperPoints is computed with the given information. \param
[in] baseGrid - base points for the HyperLines \param [in] data - contains
information about the dataset, domain, sampling and output files \param [in]
search - contains information about the considered time intervals \param [in]
mainDir - axis parallel direction, in which the search takes place has to be one
of the following: \li (1,0,0) - for grid in yz-plane, search in x-direction \li
(0,1,0) - for grid in xz-plane, search in y-direction \li (0,0,1) - for grid in
xy-plane, search in z-direction \return found recirculation points
*/
vector<RecPoint>
JobWorker::searchByBaseGrid(const vector<Vec3r>& baseGrid,
                            const DataParams&    data,
                            const SearchParams&  search,
                            const Vec3i&         mainDir,
                            const real&          prec,
                            Flow3D*              p_flow,
                            list<RecPoint>*      p_candidates) const {
  // compute the different directions in parallel
  //---> variables used for different threads (critical section)
  std::vector<RecPoint> recPoints;
  Vec3r                 lineDir =
    computeOffsetVector(data); // offset between two samples in x,y,z-direction
  lineDir = lineDir * mainDir; // offset in search direction

  FlowSampler3D    flowSampler(*p_flow);
  vector<RecPoint> threadRecPoints;
  list<RecPoint>   threadCandidates;

  for (size_t i = 0; i < baseGrid.size(); i++) {
    // get the points, which are used for the hyperlines
    Vec3r     pA = baseGrid[i];
    Vec3r     pB = pA + lineDir;
    HyperLine hl(pA, pB, &flowSampler);

    utils::printProgress(i + 1,
                         baseGrid.size(),
                         "recirculation points found: " +
                           to_string(recPoints.size()));

    // pA and pB are consecutively propagated, until the
    // whole domain in params.lineDir has been processed
    while (pB[0] <= data.dMax[0] && pB[1] <= data.dMax[1] &&
           pB[2] <= data.dMax[2]) {
      // check if the hyperline is inside the domain, if not go on
      if (!p_flow->isInside(pA) || !p_flow->isInside(pB)) {
        // move the HyperLine to the next position
        pA = pB;
        pB += lineDir;
        hl = HyperLine(hl.getHyperPointB(), pB);
        continue;
      }

      // search the HyperLine for Recirculation Points
      vector<RecPoint> tempRecPoints;
      tempRecPoints = hl.getRecirculationPoints(search.t0_min,
                                                search.t0_max,
                                                search.tau_min,
                                                search.tau_max,
                                                search.dt,
                                                prec,
                                                false);
      for (auto rcp : tempRecPoints)
        recPoints.push_back(rcp);

      // move the HyperLine to the next position
      pB += lineDir;
      hl = HyperLine(hl.getHyperPointB(), pB);
    }
  }

  std::cout << "searched " << baseGrid.size() << " items. found "
            << recPoints.size() << " RecPoints. \n";
  if (p_candidates)
    for (auto rcp : recPoints)
      (*p_candidates).push_back(rcp);

  return recPoints;
}
}
