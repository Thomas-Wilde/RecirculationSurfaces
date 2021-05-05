#include "amiradataset.hh"
#include "doublegyre3D.hh"
#include "flow.hh"
#include "globals.hh"
#include "jobparams.hh"
#include "jobworker.hh"
#include "types.hh"

#include <filesystem>
#include <iostream>
#include <memory>
#include <set>
#include <stdio.h>

using namespace std;
using namespace RS;
namespace fs = std::filesystem;

//--------------------------------------------------------------------------//
void
sampleRecirculationPointToPathline(const RecPoint& rcp,
                                   const int&      samples,
                                   FlowSampler3D&  sampler) {
  // show the rcp and some points on the pathline
  Pathline3D pathline = sampler.sampleFlow(rcp.pos, rcp.t0, rcp.tau);
  cout << "Sample pathline of Recirculation Point.\n";
  cout << "Recirculation Point at: " << rcp.toString() << "\n";

  RS::real dt = rcp.tau / static_cast<float>(samples - 1);
  for (int i = 0; i < samples; i++)
    cout << "pos: " << pathline.eval_position_at(rcp.t0 + dt * i)
         << " tau=" << dt * i << "\n";

  auto p0   = pathline.eval_position_at(rcp.t0);
  auto ptau = pathline.eval_position_at(rcp.t0 + dt * (samples - 1));
  cout << "\n\ndistance (𝕡0, 𝕡τ): " << (p0 - ptau).norm2() << endl;
}

//--------------------------------------------------------------------------//
void
computeDoubleGyre3DSingleLine(void) {
  /*********************************************/
  //-- 3D DoubleGyre place single HyperLine  --//
  /*********************************************/
  shared_ptr<Flow3D> flow = make_shared<DoubleGyre3D>();
  FlowSampler3D      sampler(*flow);

  //--- Set 2
  Vec3r pA(1.6708, 0.538, 0.475);
  Vec3r pB(1.6708, 0.538, 0.525);

  // check if the hyperline is inside the domain, if not go on
  if (!flow->isInside(pA) || !flow->isInside(pB))
    return;

  // search the HyperLine for Recirculation Points
  HyperLine        hl(pA, pB, &sampler);
  vector<RecPoint> recPoints;
  recPoints =
    hl.getRecirculationPoints(0.0, 10.0, 0.0, 10.0, 0.2, Globals::SEARCHPREC);

  for (auto rcp : recPoints)
    sampleRecirculationPointToPathline(rcp, 20, sampler);
}

//--------------------------------------------------------------------------//
void
computeDoubleGyre3D(void) {
  /*********************************************/
  //---- 3D DoubleGyre full domain tau = 10 ---//
  /*********************************************/
  Vec3r dMin(0.01, 0.01, 0.01);
  Vec3r dMax(1.99, 0.99, 0.99);
  Vec3i sampling(201, 101, 101);

  RS::real prec = Globals::SEARCHPREC;

  DataParams         data(dMin, dMax, sampling, "DoubleGyre3D", "");
  SearchParams       search(0., 10., 0., 15., 0.2, prec * 100., prec / 100.);
  shared_ptr<Flow3D> flow = make_shared<DoubleGyre3D>();

  JobWorker jw;

  // start search: results are saved automatically
  data.saveFilename = "DoubleGyre_xDir.ply";
  jw.searchDomainByBaseGrid(data, search, Vec3i(1, 0, 0), flow.get());

  data.saveFilename = "DoubleGyre_yDir.ply";
  jw.searchDomainByBaseGrid(data, search, Vec3i(0, 1, 0), flow.get());

  data.saveFilename = "DoubleGyre_zDir.ply";
  jw.searchDomainByBaseGrid(data, search, Vec3i(0, 0, 1), flow.get());
}

//--------------------------------------------------------------------------//
void
computeSquaredCylinderSingleLine(void) {
  // ---------------------------------------------------------------------- //
  // --- load the data set
  string      path = "../../data/SquareCylinderHighResTime";
  set<string> file_set; // we need this container for sorting
  for (const auto& entry : fs::directory_iterator(path))
    if (entry.is_regular_file())
      file_set.insert(entry.path());
  // ---
  std::vector<std::string> file_vec; // the amira flow needs a vector
  for (const auto& entry : file_set) {
    file_vec.push_back(entry);
    if (60 == file_vec.size()) // control how many files get loaded
      break;                   // for the example we load data up to t=7.0
  }

  // ---------------------------------------------------------------------- //
  // --- create a flow based on the amira files
  // Vec2r time_range(0.08, 40.56);  // you get this info from the filenames
  Vec2r time_range(0.0, 0.08 * file_vec.size());
  Vec3r base_pnt(0.5, -0.5, 0.0); // this is a point directly behind the
                                  // cylinder  we only have to search
                                  // behind this point
  shared_ptr<Flow3D> flow = make_shared<AmiraDataSet>(file_vec, time_range);
  FlowSampler3D      sampler(*flow);
  cout << "time range: [" << time_range << "]" << endl;

  // ----------------------------------------------------------------------
  //
  //--- We know a recirculation point with these coordinates
  // (x,y,z,t,τ) = (1.42 0.46 2.01549 1.22299 3.84389)
  Vec3r pA(1.42, 0.46, 2.00);
  Vec3r pB(1.42, 0.46, 2.02);

  //--- check if the hyperline is inside the domain, if not go on
  if (!flow->isInside(pA) || !flow->isInside(pB))
    return;

  // search the HyperLine for Recirculation Points
  HyperLine        hl(pA, pB, &sampler);
  vector<RecPoint> recPoints;
  recPoints = hl.getRecirculationPoints(0.0, 5.0, 0.0, 5.0, 0.2, 0.0005);
  for (auto rcp : recPoints)
    sampleRecirculationPointToPathline(rcp, 20, sampler);
}

//--------------------------------------------------------------------------//
int
main(int /*argc*/, char** /*argv[]*/) {
  cout << setprecision(8) << fixed << showpos << "\n";
  // computeDoubleGyre3D();
  // computeDoubleGyre3DSingleLine();
  computeSquaredCylinderSingleLine();
  return 0;
}
