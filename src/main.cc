#include <iostream>
#include <memory>
#include <stdio.h>

#include "doublegyre3D.hh"
#include "flow.hh"
#include "globals.hh"
#include "jobparams.hh"
#include "jobworker.hh"
#include "types.hh"

using namespace std;
using namespace RS;

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
  for (int i = 0; i < 20; i++)
    cout << "pos: " << pathline.eval_position_at(rcp.t0 + dt * i)
         << " tau=" << dt * i << "\n";
}

//-----------------------------------------------------------------------------------------------//
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

int
main(int /*argc*/, char** /*argv[]*/) {
  cout << setprecision(8) << fixed << showpos << "\n";
  // computeDoubleGyre3D();
  computeDoubleGyre3DSingleLine();
  return 0;
}
