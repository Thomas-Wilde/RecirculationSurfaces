#ifndef _H_JOBWORKER_
#define _H_JOBWORKER_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>

#include "flow.hh"
#include "hyperline.hh"
#include "iorecpointfile.hh"
#include "jobparams.hh"
#include "recpoint.hh"

// #include "doublegyre.h"
// #include "doublegyre3D.h"
// #include "flowgenerator.h"
// #include "hypercube.h"
// #include "icsquad.h"
// #include "steadydoublegyre3d.h"
// #include "vectorfieldcuboid.h"

namespace RS {
class JobWorker {
public:
  JobWorker(void);
  ~JobWorker(void);

  void searchDomainByBaseGrid(
    const DataParams&      data,
    const SearchParams&    search,
    const Vec3i&           mainDir,
    Flow3D*                p_flow      = nullptr,
    std::vector<RecPoint>* p_recPoints = nullptr) const;

  bool generateBaseGrid(std::vector<Vec3r>* p_baseGrid,
                        const DataParams&   data,
                        const Vec3i&        searchDir) const;

  Vec3r computeOffsetVector(const DataParams& data) const;

  std::vector<RecPoint> searchByBaseGrid(
    const std::vector<Vec3r>& baseGrid,
    const DataParams&         data,
    const SearchParams&       search,
    const Vec3i&              mainDir,
    const real&               prec,
    Flow3D*                   p_flow       = nullptr,
    std::list<RecPoint>*      p_candidates = nullptr) const;

  //   void start(void);

  //   void searchDomainByLines(const DataParams& data, const SearchParams&
  //   search);

  //   bool preparePlyFile(std::fstream*      p_fileStream,
  //                       const std::string& filename,
  //                       const int&         pointCount = 0) const;
  //   bool writeRecPointsToStream(std::fstream*                p_fileStream,
  //                               const std::vector<RecPoint>& recPoints,
  //                               bool onlyHits = true) const;
  //   void generateSamples(std::vector<Vec3r>* p_sampleVec,
  //                        const DataParams&   data) const;

  //   void computeClosedStreamline(flow::Flow3D&     flow,
  //                                const DataParams& data,
  //                                const Vec3i&      sampling,
  //                                rs::real          x_pos);
  //   void computeClosedStreamline2(flow::Flow3D&     flow,
  //                                 const DataParams& data,
  //                                 const Vec3i&      sampling,
  //                                 rs::real          x_pos);

  //   std::vector<Vec4r> computeClosedStreamline3(flow::Flow3D&     flow,
  //                                               const DataParams& data,
  //                                               const Vec3i&      dir,
  //                                               const rs::real&   pos,
  //                                               const rs::real&   tau_min,
  //                                               const rs::real&   tau_max,
  //                                               const rs::real&   dt) const;

  //   ///--->
  //   //--- BaseGrid ---//
  // public:
  //   //--- adative search ---//
  //   void                  searchDomainAdaptive(const DataParams&      data,
  //                                              const SearchParams&    search,
  //                                              const Vec3i& mainDir,
  //                                              std::vector<RecPoint>*
  //                                              p_recPoints = nullptr) const;
  //   std::vector<RecPoint> refinePointsOnLine(const DataParams&      data,
  //                                            const SearchParams&    search,
  //                                            std::vector<RecPoint>&
  //                                            candidates, const Vec3i&
  //                                            mainDir, const real& offset)
  //                                            const;

  //   //--- refinment of given recirculation points ---//
  //   std::vector<RecPoint> refinePointsSpace(
  //     const DataParams&            data,
  //     const SearchParams&          search,
  //     const std::vector<RecPoint>& recPoints,
  //     flow::Flow3D*                p_flow = nullptr) const;
  //   std::vector<Vec3r> generateRefineLines(const Vec3r&      rcpPos,
  //                                          const DataParams& data) const;

  // private:
  //   std::vector<RecPoint> searchByBaseGrid(
  //     const std::vector<Vec3r>& baseGrid,
  //     const DataParams&         data,
  //     const SearchParams&       search,
  //     const Vec3i&              mainDir,
  //     const real&               prec,
  //     flow::Flow3D*             p_flow       = nullptr,
  //     std::list<RecPoint>*      p_candidates = nullptr) const;

  //   void superviseSearchBaseGrid(const int&                threadID,
  //                                const int&                numThreads,
  //                                const std::vector<int>&   processedItems,
  //                                const std::vector<int>&   foundPoints,
  //                                const std::vector<int>&   currentItem,
  //                                std::vector<int>*         p_repeatItem,
  //                                std::vector<bool>*        p_finished,
  //                                std::vector<int>*         p_endProcess,
  //                                const std::string&        filename,
  //                                const std::vector<Vec3r>& baseGrid,
  //                                const int&                maxRepeat,
  //                                bool* p_init_supervise) const;

  //   void superviseRefinePointsSpace(const int&              threadID,
  //                                   const int&              numThreads,
  //                                   const std::vector<int>& processedItems,
  //                                   const std::vector<int>& foundPoints,
  //                                   const std::vector<int>& currentItem,
  //                                   std::vector<int>*       p_repeatItem,
  //                                   std::vector<bool>*      p_finished,
  //                                   std::vector<int>*       p_endProcess,
  //                                   const std::string&      filename,
  //                                   const int&              itemCount,
  //                                   const int&              maxRepeat,
  //                                   bool* p_init_superviser) const;

  //   bool prepareOutputFiles(
  //     const std::string& filename, // filename without extension
  //     std::fstream* p_pointFile, // file to save the compute recirculation
  //     points std::fstream* p_logFile) const;

  // private:
  //   void sendMessage(QString message) const;

  //   std::vector<RecPoint> searchAtSamples(
  //     const std::vector<Vec3r>& samples,
  //     const DataParams&         data,
  //     const SearchParams&       params,
  //     const Vec3r&              searchDir,
  //     real*                     p_integrateTime   = nullptr,
  //     real*                     p_trilinearTime   = nullptr,
  //     std::map<int, int>*       p_recDepthMapHits = nullptr,
  //     std::map<int, int>*       p_recDepthMapMiss = nullptr,
  //     unsigned int*             p_integrateCalls  = nullptr);

  //   bool savePositions(const std::vector<Vec3r>& posVec,
  //                      const std::string&        filename) const;

  //

  //   int m_sampleCount;
};
} // namespace rs
#endif //_H_JOBWORKER_
