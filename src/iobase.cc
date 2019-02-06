#include "iobase.hh"
#include <cstring>

namespace RS {

//--------------------------------------------------------------------------//
IOBase::IOBase(void) {}

//--------------------------------------------------------------------------//
IOBase::~IOBase(void) {}

//--------------------------------------------------------------------------//
const std::string&
IOBase::getInfo(void) const {
  return m_info;
}

//--------------------------------------------------------------------------//
bool
IOBase::openFile(const char* p_filename, const char* p_mode) {
  // check if a file is already opened
  if (mp_file != nullptr) {
    m_info = std::string("Tried to open more than one file at a time.");
    return false;
  }

  // open the stream
  mp_file = fopen(p_filename, p_mode);
  if (!mp_file) { // check if something went wrong, e.g. file does not exist
    m_info = std::string("Could not open ");
    m_info.append(p_filename);
    return false;
  }

  m_info = std::string("Opened file ");
  m_info.append(p_filename);
  return true;
}

//--------------------------------------------------------------------------//
bool
IOBase::closeFile(void) {
  if (mp_file == nullptr) {
    m_info = std::string("Tried to close file but none is opened.");
    return false;
  }

  fclose(mp_file);
  m_info  = std::string("Closed file. \n");
  mp_file = nullptr;
  return true;
}

//--------------------------------------------------------------------------//
bool
IOBase::exitWithFailure(void) {
  // closing the file leads to overwriting the info string so we temporarily
  // save it
  std::string tempInfo = m_info;
  closeFile();
  m_info.append(tempInfo);
  return false;
}

//--------------------------------------------------------------------------//
bool
IOBase::exitWithFailure(const std::string& errorMessage) {
  closeFile();
  m_info = errorMessage;
  return false;
}

//--------------------------------------------------------------------------//
bool
IOBase::exitWithFailure(const char* p_info, const char* p_add) {
  std::string tempString = p_info;
  tempString.append(p_add);
  return exitWithFailure(tempString);
}

//--------------------------------------------------------------------------//
/** Find a string in the given buffer and return a pointer
to the contents directly behind the SearchString.
If not found, return the buffer. A subsequent sscanf()
will fail then, but at least we return a decent pointer. */
bool
IOBase::findAndJump(const char*  p_buffer,
                    const char*  p_searchString,
                    const char** p_pos) const {
  const char* FoundLoc = strstr(p_buffer, p_searchString);
  if (FoundLoc) {
    if (p_pos != nullptr)
      *p_pos = FoundLoc + strlen(p_searchString);
    return true;
  }
  return false;
}

//--------------------------------------------------------------------------//
/**\return \li true, if the current line matches 'p_compare'
           \li false, if the string does not match*/
bool
IOBase::checkContent(const char* p_buffer, const char* p_compare) const {
  const char* p_pos = &p_buffer[0];
  if (findAndJump(p_buffer,
                  p_compare,
                  &p_pos)) // check if the string is somewhere in the buffer
    return (p_pos - strlen(p_compare) ==
            &p_buffer[0]); // return true, if it is at the beginning

  return false;
}
}
