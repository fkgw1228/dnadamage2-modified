#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(dnadamage2_modified)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#-------------------------------------------------------------------------------
# Eigen3 and nlohmann_json
#
find_package(Eigen3 3.3 REQUIRED)
find_package(nlohmann_json REQUIRED)

#----------------------------------------------------------------------------
# Add the executable, use our local headers, and link it to the Geant4 libraries
#
add_executable(dnadamage2_modified dnadamage2_modified.cc ${sources} ${headers})
target_include_directories(dnadamage2_modified PRIVATE include)
target_link_libraries(dnadamage2_modified
  PRIVATE
    ${Geant4_LIBRARIES}
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB dnadamage2_SCRIPTS
  ${PROJECT_SOURCE_DIR}/*.in
  ${PROJECT_SOURCE_DIR}/*.mac
  ${PROJECT_SOURCE_DIR}/*.txt
  ${PROJECT_SOURCE_DIR}/*.C
  )
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/analysis DESTINATION ${CMAKE_BINARY_DIR})

foreach(_script ${dnadamage2_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dnadamage2_modified DESTINATION bin )
