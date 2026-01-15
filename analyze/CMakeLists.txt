cmake_minimum_required(VERSION 3.15)

project(dnadamage2-modified-plotter)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT library
find_package(ROOT REQUIRED COMPONENTS
    RIO
    Hist
    Graf
    Graf3d
)

# Eigen3 library
find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)

# Include directories
include_directories(
    ${EIGEN3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add executables
add_executable(plot_damage
    PlotDamageLocation.cc
)
add_executable(analyze_damage
    AnalyzeDamage.cc
    DamageAnalyzer.cc
)
add_executable(report_damage
    DamageAnalyzer.cc
    SBRatioReporter.cc
)

file(COPY PlotSsbDsb.C DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Link libraries
target_link_libraries(plot_damage
    PRIVATE
        ROOT::RIO
        ROOT::Hist
        ROOT::Graf
        ROOT::Graf3d
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
)
