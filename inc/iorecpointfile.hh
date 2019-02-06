#ifndef _H_IORECPOINTFILE_
#define _H_IORECPOINTFILE_

#include "iobase.hh"
#include "recpoint.hh"
#include <vector>

namespace RS {
class IORecPointFile : public IOBase {
public:
  IORecPointFile(void);
  virtual ~IORecPointFile(void);

  /**Load recircualtion points from file and append them to p_recPoints.
     \param p_recPoints - the points read from file are \b appended on this
     vector*/
  bool loadFile(const std::string&     filename,
                std::vector<RecPoint>* p_recPoints);
  /**Load recircualtion points and triangulation. The data in p_recPoints and
    p_triangles is \b replaced. \param p_recPoints - the points read from file
    are stored here, previous data is deleted \param p_triangles - the indexes
    of the points forming the triangles are stored here, previous data is
    deleted*/
  bool loadFile(const std::string&     filename,
                std::vector<RecPoint>* p_recPoints,
                std::vector<int>*      p_triangles);
  bool saveFile(const std::string&           filename,
                const std::vector<RecPoint>& recPoints,
                std::string                  author = "---");

  bool loadMeshlabFile(const std::string&  filename,
                       std::vector<Vec3r>* p_points,
                       std::vector<int>*   p_triangles);

protected:
  bool saveFile(const std::string& /*filename*/) override { return false; }
  bool loadFile(const std::string& /*filename*/) override { return false; };
  int  loadRecPoints(std::vector<RecPoint>* p_recPoints);
  int  loadPoints(std::vector<Vec3r>* p_points);
  int  loadTriangles(std::vector<int>* p_triangles);

  bool         loadHeader(bool skip_last_lines = true);
  unsigned int m_rcpCount;
  unsigned int m_triangles_count = 0;
  char         m_readBuffer[2048];
};
}

#endif //_H_IORECPOINTFILE_
