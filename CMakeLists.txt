# At LEAST 2.8 but newer is better
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(RecirculationSurfaces VERSION 0.1 LANGUAGES CXX)

# Set the path where cmake modules are stored
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
message("Module path: " ${CMAKE_MODULE_PATH})

#export compilation database to use with auto complete
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# search for libraries that are needed: Flann and LZ4
find_package(LZ4 REQUIRED)

if (LZ4_FOUND)
  message("LZ4 include directory: " STATUS ${LZ4_INCLUDE_DIR})
else(LZ4_FOUND)
  message("ERROR: LZ4 could not be found.")
endif()

find_package(Flann REQUIRED)
if (FLANN_FOUND)
  message("Flann include directory:" STATUS ${FLANN_INCLUDE_DIR})
  message("Flann libraries:" STATUS ${FLANN_LIBRARIES})
else(FLANN_FOUND)
  message("ERROR: Flann could not be found.")
endif()

# Define an executable
add_executable(recsurface
               src/main.cc
               src/globals.cc
               src/jobworker.cc
               src/amiradataset.cc
               src/doublegyre3D.cc
               src/hyperline.cc
               src/hyperpoint.cc
               src/critextractor.cc
               src/vectorcuboid.cc
               src/iorecpointfile.cc
               src/iobase.cc
               vclibs/math/rk43.cc
               vclibs/math/ode.cc)

target_include_directories(
  recsurface
    PUBLIC ${FLANN_INCLUDE_DIRS}
    PUBLIC ${LZ4_INCLUDE_DIR}
    PRIVATE ./inc
    PRIVATE ./)

# dependencies for libraries
target_link_libraries(
  recsurface
    PUBLIC ${FLANN_LIBRARIES}
    PUBLIC ${LZ4_LIBRARY}
    PRIVATE stdc++fs)

# additional compiler flags for more warnings
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -llz4")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -llz4")

