#ifndef _H_TYPES_
#define _H_TYPES_

//#include "vclibs/math/Mat.hh"
#include "vclibs/math/VecN.hh"
#include "vclibs/math/rk43.hh"

namespace RS {

typedef double       real;
typedef unsigned int unInt;

typedef VC::math::VecN<real, 2> Vec2r;
typedef VC::math::VecN<real, 3> Vec3r;
typedef VC::math::VecN<real, 4> Vec4r;
typedef VC::math::VecN<real, 5> Vec5r;
typedef VC::math::VecN<real, 6> Vec6r;
typedef VC::math::VecN<int, 2>  Vec2i;
typedef VC::math::VecN<int, 3>  Vec3i;
typedef VC::math::VecN<int, 4>  Vec4i;
typedef VC::math::VecN<int, 5>  Vec5i;
typedef VC::math::VecN<int, 6>  Vec6i;
// typedef VC::math::Mat<real, 2, 2> Mat2r;

typedef VC::math::ode::Solution<real, Vec3r> Pathline3D;
typedef VC::math::ode::Solution<real, Vec2r> Pathline2D;

struct Range {
  Range(real min_val, real max_val)
  : min(min_val)
  , max(max_val) {}

  Range(void)
  : min(-std::numeric_limits<double>::max())
  , max(std::numeric_limits<double>::max()) {}

  real min = -std::numeric_limits<double>::max();
  real max = std::numeric_limits<double>::max();

  operator std::string() const {
    return std::string("[") + std::to_string(min) + std::string(", ") +
           std::to_string(max) + std::string("]");
  }
};

struct Domain {
  Domain(const Range& r_x, const Range& r_y, const Range& r_z, const Range& r_t)
  : range_x(r_x)
  , range_y(r_y)
  , range_z(r_z)
  , range_t(r_t) {}

  Domain(const real& x_min,
         const real& x_max,
         const real& y_min,
         const real& y_max,
         const real& z_min,
         const real& z_max,
         const real& t_min,
         const real& t_max)
  : range_x(x_min, x_max)
  , range_y(y_min, y_max)
  , range_z(z_min, z_max)
  , range_t(t_min, t_max) {}

  Range range_x;
  Range range_y;
  Range range_z;
  Range range_t;
};

template<typename T>
struct CmpVec {
  CmpVec<T>(void){};

  /**Compares elements in a and b componentwise lexicographically. The order in
     which the dimensions are compared can be controlled by the permutation
     vector. If e.g. a and b are three dimensional vectors and \c p_perm_vec ==
     nullptr, then the components of a and b are compared in xyz-order, If \c
     p_perm_vec is given with [1,0,2], then the vectors are compared in
     yxz-order. \throws out-of-range exception if the permutation vector
     contains indices, which are too big or if the vector itself smaller than
     the length of a. \param p_perm_vec - permutation vector with 0-based
     indexes \return \li true if a < b (sorted lexicographically)
     \li false if a >= b */
  bool operator()(T const&                         a,
                  T const&                         b,
                  const std::vector<unsigned int>* p_perm_vec = nullptr) const {
    real dif   = 0.0;
    real sum   = 0.0;
    bool equal = false;
    real epsi  = std::numeric_limits<double>::epsilon();
    real prec  = 10.0;

    // check for equality until the first difference
    unsigned int i     = 0;
    unsigned int index = 0;
    for (i = 0; i < a.size(); ++i) {
      index = p_perm_vec ? (*p_perm_vec)[i] : i;
      if (index >= a.size())
        throw std::out_of_range(
          "Index is out of bounds. Check permutation vector!");
      dif = a[index] - b[index];
      sum = a[index] + b[index];
      equal =
        (std::abs(dif) < prec * epsi * std::abs(sum)) || std::abs(dif) < epsi;
      if (!equal)
        break; // search only to the first difference
    }
    if (i == a.size())
      return false; // this implies (a==b)
    return dif < 0.0;
  }

  bool operator()(const T*                         p_a,
                  const T*                         p_b,
                  const std::vector<unsigned int>* p_perm_vec = nullptr) const {
    return operator()(*p_a, *p_b, p_perm_vec);
  }
};

struct CmpVec2i : CmpVec<Vec2i> {};
struct CmpVec3i : CmpVec<Vec3i> {};
struct CmpVec2r : CmpVec<Vec2r> {};
struct CmpVec3r : CmpVec<Vec3r> {};
struct CmpVec4r : CmpVec<Vec4r> {};

} // namespace RS
#endif //_H_TYPES_
