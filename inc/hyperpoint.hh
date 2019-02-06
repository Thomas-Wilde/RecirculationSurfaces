#ifndef _H_HYPERPOINT_
#define _H_HYPERPOINT_

#include <map>
#include <memory>

#include "flowsampler.hh"
#include "types.hh"

namespace RS {

class HyperPoint {
public:
  HyperPoint(const Vec3r& pos, RS::FlowSampler3D* p_flowSampler);
  HyperPoint(const HyperPoint& point);
  ~HyperPoint(void);
  const std::shared_ptr<FlowMap3D> getFlowMap(const real& t0, const real& tau);
  void                             addFlowMap(const real&                t0,
                                              const real&                tau,
                                              std::shared_ptr<FlowMap3D> flow_map);

  HyperPoint&    operator=(const HyperPoint& point);
  FlowSampler3D* getFlowSampler(void) const;
  void           setFlowSampler(FlowSampler3D* p_flowSampler);

  const Vec3r& getPos(void) const;

private:
  Vec3r          m_pos;
  FlowSampler3D* mp_flowSampler;

  std::map<std::pair<real, real>, std::shared_ptr<FlowMap3D>> m_flowMaps;
};
} // namespace rs

#endif //_H_HYPERPOINT_
