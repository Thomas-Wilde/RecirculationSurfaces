#include "iorecpointfile.hh"

namespace RS {
//-----------------------------------------------------------------------------------------------//
IORecPointFile::IORecPointFile(void)
: m_rcpCount(0) {}

//-----------------------------------------------------------------------------------------------//
IORecPointFile::~IORecPointFile(void) {}

//-----------------------------------------------------------------------------------------------//
bool
IORecPointFile::loadFile(const std::string&     filename,
                         std::vector<RecPoint>* p_recPoints,
                         std::vector<int>*      p_triangles) {
  assert(p_recPoints != nullptr);
  //--- try to open the file
  if (!openFile(&filename[0], "r"))
    return false;
  //--- analyze header and jump to data section
  if (!loadHeader())
    return exitWithFailure();
  //--- load the points
  p_recPoints->clear();
  unsigned int c = loadRecPoints(p_recPoints);
  // check for amount of read data
  if (m_rcpCount != c)
    return exitWithFailure("Something went wrong while reading the point data "
                           "section.\nPremature end of file?\n");
  //--- load the triangles
  c = loadTriangles(p_triangles);
  // check for amount of read data
  if (m_triangles_count != c)
    return exitWithFailure("Something went wrong while reading the triangle "
                           "data section.\nPremature end of file?\n");

  closeFile();
  m_info = "Succesfully loaded " + std::to_string(c) + " recirculation points";
  return true;
}

//-----------------------------------------------------------------------------------------------//
bool
IORecPointFile::loadFile(const std::string&     filename,
                         std::vector<RecPoint>* p_recPoints) {
  assert(p_recPoints != nullptr);
  //--- try to open the file
  if (!openFile(&filename[0], "r"))
    return false;
  //--- analyze header and jump to data section
  if (!loadHeader())
    return exitWithFailure();
  //--- load the points
  unsigned int c = loadRecPoints(p_recPoints);
  // check for amount of read data
  if (m_rcpCount != c)
    return exitWithFailure("Something went wrong while reading the data "
                           "section.\nPremature end of file?\n");
  closeFile();
  m_info = "Succesfully loaded " + std::to_string(c) + " recirculation points";
  return true;
}

//-----------------------------------------------------------------------------------------------//
int
IORecPointFile::loadRecPoints(std::vector<RecPoint>* p_recPoints) {
  // load the RecPoints
  unsigned int c = 0;
  real         tempArray[6];
  real*        p       = &tempArray[0];
  unsigned int tempSID = 0;
  p_recPoints->reserve(p_recPoints->size() + m_rcpCount);
  while (c < m_rcpCount) {
    fgets(m_readBuffer, 2047, mp_file);
    RecPoint rcp;
    sscanf(&m_readBuffer[0],
           "%lf %lf %lf %lf %lf %lf %u",
           &p[0],
           &p[1],
           &p[2],
           &p[3],
           &p[4],
           &p[5],
           &tempSID);
    rcp.pos[0] = tempArray[0];
    rcp.pos[1] = tempArray[1];
    rcp.pos[2] = tempArray[2];
    rcp.t0     = tempArray[3];
    rcp.tau    = tempArray[4];
    rcp.dist   = tempArray[5];
    rcp.sId    = tempSID;
    p_recPoints->push_back(rcp);
    ++c;
  }
  return c;
}

//-----------------------------------------------------------------------------------------------//
int
IORecPointFile::loadTriangles(std::vector<int>* p_triangles) {
  // load the RecPoints
  unsigned int c = 0;
  int          temp_array[4];
  int*         p = &temp_array[0];
  p_triangles->clear();
  p_triangles->reserve(m_triangles_count * 3);
  while (c < m_triangles_count) {
    fgets(m_readBuffer, 2047, mp_file);
    // sscanf_s(&m_readBuffer[0], "%d %d %d %d",
    //          &p[0], &p[1], &p[2], &p[3] );
    sscanf(&m_readBuffer[0], "%d %d %d %d", &p[0], &p[1], &p[2], &p[3]);
    (*p_triangles).push_back(temp_array[1]);
    (*p_triangles).push_back(temp_array[2]);
    (*p_triangles).push_back(temp_array[3]);
    ++c;
  }
  return c;
}

//-----------------------------------------------------------------------------------------------//
bool
IORecPointFile::loadHeader(bool skip_last_lines) {
  fread(m_readBuffer, sizeof(char), 2047, mp_file);
  m_readBuffer[2047] = '\0';

  //--- check for correct header
  const char* pos = &m_readBuffer[0];
  if (!checkContent(pos, "ply")) {
    m_info = "First line should be 'ply'.";
    return false;
  }
  findAndJump(pos, "\n", &pos);

  if (!checkContent(pos, "format ascii 1.0")) {
    m_info = "Second line should be 'format ascii 1.0'";
    return false;
  }
  findAndJump(pos, "\n", &pos);

  if (!checkContent(pos, "comment content Recirculation Points") &&
      !checkContent(pos, "comment VCGLIB generated")) {
    m_info = "Third line should be 'comment content Recirculation Points' or "
             "'comment VCGLIB generated'";
    return false;
  }
  findAndJump(pos, "\n", &pos);

  // read point count
  const char* p_pos = nullptr;
  if (!findAndJump(m_readBuffer, "element vertex ", &p_pos)) {
    m_info = "Could not read element vertex count";
    return false;
  }
  sscanf(p_pos, "%u", &m_rcpCount);
  if (0 == m_rcpCount) {
    m_info = "Element vertex count is set to 0.";
    return false;
  }
  findAndJump(pos, "\n", &pos);

  // check if triangles are given
  if (findAndJump(m_readBuffer, "element face ", &p_pos)) {
    // read the amount of triangle data
    sscanf(p_pos, "%u", &m_triangles_count);
    if (0 == m_triangles_count) {
      m_info = "Triangles were expected, but count is 0.";
      return false;
    }
    findAndJump(p_pos, "\n", &p_pos);

    // check the property section for triangles
    if (!checkContent(p_pos, "property list uchar int vertex_indices")) {
      m_info = "The line after face count should be 'property list uchar int "
               "vertex_indices'";
      return false;
    }
  }

  // Find the beginning of the data section
  if (!findAndJump(pos, "end_header", &pos)) {
    m_info = "Could not find 'end_header'.";
    return false;
  }
  findAndJump(pos, "\n", &pos);

  // Set the file pointer to the beginning of "# Data section follows"
  long idxStartData = static_cast<long>(pos - &m_readBuffer[0]);
  fseek(mp_file, idxStartData, SEEK_SET);
  // consume this line, which is "end_header"
  if (skip_last_lines) {
    fgets(m_readBuffer, 2047, mp_file);
    fgets(m_readBuffer, 2047, mp_file);
  }
  return true;
}

//-----------------------------------------------------------------------------------------------//
bool
IORecPointFile::saveFile(const std::string&           filename,
                         const std::vector<RecPoint>& recPoints,
                         std::string                  author) {
  //--- try to open the file
  if (!openFile(&filename[0], "w"))
    return false;

  //--- write the header
  std::string headerPart1 = "ply\nformat ascii 1.0\ncomment content "
                            "Recirculation Points\ncomment author ";
  std::string headerPart2 = author;
  std::string headerPart3 = "\nelement vertex ";
  std::string headerPart4 = std::to_string(recPoints.size());
  std::string headerPart5 =
    "\nproperty double x\nproperty double y\nproperty double z\nproperty "
    "double t0\nproperty double tau\nend_header\n";

  fputs(&headerPart1[0], mp_file);
  fputs(&headerPart2[0], mp_file);
  fputs(&headerPart3[0], mp_file);
  fputs(&headerPart4[0], mp_file);
  fputs(&headerPart5[0], mp_file);

  for (auto rcp : recPoints) {
    fprintf(mp_file, "%.8lf ", rcp.pos[0]);
    fprintf(mp_file, "%.8lf ", rcp.pos[1]);
    fprintf(mp_file, "%.8lf ", rcp.pos[2]);
    fprintf(mp_file, "%.8lf ", rcp.t0);
    fprintf(mp_file, "%.8lf\n", rcp.tau);
  }
  closeFile();
  return true;
}

//-----------------------------------------------------------------------------------------------//
bool
IORecPointFile::loadMeshlabFile(const std::string&  filename,
                                std::vector<Vec3r>* p_points,
                                std::vector<int>*   p_triangles) {
  assert(p_points != nullptr);
  //--- try to open the file
  if (!openFile(&filename[0], "r"))
    return false;
  //--- analyze header and jump to data section
  if (!loadHeader(false))
    return exitWithFailure();
  //--- load the points
  p_points->clear();
  unsigned int c = loadPoints(p_points);
  // check for amount of read data
  if (m_rcpCount != c)
    return exitWithFailure("Something went wrong while reading the point data "
                           "section.\nPremature end of file?\n");
  //--- load the triangles
  c = loadTriangles(p_triangles);
  // check for amount of read data
  if (m_triangles_count != c)
    return exitWithFailure("Something went wrong while reading the triangle "
                           "data section.\nPremature end of file?\n");

  closeFile();
  m_info = "Succesfully loaded " + std::to_string(c) + " recirculation points";
  return true;
}

//-----------------------------------------------------------------------------------------------//
int
IORecPointFile::loadPoints(std::vector<Vec3r>* p_points) {
  // load the RecPoints
  unsigned int c = 0;
  real         tempArray[3];
  real*        p = &tempArray[0];
  p_points->reserve(p_points->size() + m_rcpCount);
  while (c < m_rcpCount) {
    fgets(m_readBuffer, 2047, mp_file);
    Vec3r pnt;
    sscanf(&m_readBuffer[0], "%lf %lf %lf", &p[0], &p[1], &p[2]);
    pnt[0] = tempArray[0];
    pnt[1] = tempArray[1];
    pnt[2] = tempArray[2];
    p_points->push_back(pnt);
    ++c;
  }
  return c;
}
} // namespace RS
