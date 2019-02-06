#ifndef _H_UTILS_
#define _H_UTILS_

#include "types.hh"

namespace RS {
namespace utils {
/**Helper function to supress warnings due to unused paramters*/
template<class... Types>
void
unusedArgs(const Types... /*args*/) {}

/**Check if the components of two vectors have the same sign.
   \return true, if all component pairs have the same sign.
   false, if at least one component pair has different signs*/
template<typename T, unsigned int n>
bool
checkSigns(const VC::math::VecN<T, n>& vec0, const VC::math::VecN<T, n>& vec1) {
  for (unsigned int i = 0; i < n; i++)
    if (vec0[i] * vec1[i] < 0.0)
      return false;
  return true;
}

//-----------------------------------------------------------------------------------------------//
/**Print a progress bar to the console. The progress bar is scaled to 20 steps.
   A message can be appended.*/

static inline void
printProgress(const int& step, const int& stepMax, const std::string msg) {
  int   scaledMax = 20;
  float scale     = static_cast<float>(scaledMax) / static_cast<float>(stepMax);
  int   stepScale =
    std::min(scaledMax - 1, static_cast<int>(static_cast<float>(step) * scale));

  std::cout << "[";
  int i = 0;
  for (i = 0; i <= stepScale; i++)
    std::cout << "x";
  for (; i < scaledMax; i++)
    std::cout << " ";
  std::cout << "] " << std::to_string(step) << "/" << std::to_string(stepMax)
            << " " << msg;
  std::cout.flush();
  std::cout << "\r"; // jump to start of the line
}
static inline void
printProgress(const int& step, const int& stepMax) {
  printProgress(step, stepMax, "");
}
} // namespace utils
} // namespace FT
#endif
