#ifndef _H_IOBASE_
#define _H_IOBASE_

#include "stdio.h"
#include <sstream>
#include <string>

namespace RS {
class IOBase {
public:
  IOBase(void);
  virtual ~IOBase(void);

  virtual bool loadFile(const std::string& filename) = 0;
  virtual bool saveFile(const std::string& filename) = 0;

  const std::string& getInfo(void) const;

  bool findAndJump(const char*  p_buffer,
                   const char*  p_searchString,
                   const char** p_pos) const;
  bool checkContent(const char* p_buffer, const char* p_compare) const;

protected:
  /**\param p_mode \li "r" - read
                   \li "w" - write
                   \li "a" - append
                   \li "rb" - read binary
                   \li "wb" - write binary
                   \li "ab" - append binary */
  bool openFile(const char* p_filename, const char* p_mode);
  bool closeFile(void);

  bool exitWithFailure(void);
  bool exitWithFailure(const std::string& errorMessage);
  bool exitWithFailure(const char* p_info, const char* p_add);

  std::string m_info  = "";
  FILE*       mp_file = nullptr;
};
} // namespace RS

#endif //_H_IOBASE_
