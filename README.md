README: [Japanese](README.ja.md)

# dnadamage2-modified

* [Overview](#overview)
* [Requirements](#requirements)
* [Installation](#installation)
  * [Installation of mandatory packages](#installation-of-mandatory-packages)
  * [Installation of this project](#installation-of-this-project)
* [Usage](#usage)

## Overview
This Geant4 application is based on the `dnadamage2` example from Geant4-DNA. It simulates DNA damage caused by ionizing radiation, covering both direct and indirect effects. It also includes analysis scripts to compute DNA damage results.

## Requirements
- Geant4 version 11.3.0 to 11.3.2
- CMake 3.16 or higher
- C++17-compatible compiler
- Python 3.8 or higher (for analysis code)

## Installation
### Installation of mandatory packages

Linux (Debian/Ubuntu):
```
sudo apt update
sudo apt install -y libeigne3-dev nlohmann-json3-dev
```
Linux (build from source)
```
# A install directory
mkdir ~/local
# Eigen3
git clone https://gitlab.com/libeigen/eigen.git
cd eigen && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/local
make install
# nlohmann_json
git clone https://github.com/nlohmann/json.git
cd json && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/local -DJSON_BuildTests=OFF
make install
# check installation
ls ~/local/include         # eigen3/  nlohmann/ will apear if they are installed
```
Mac (Homebrew):
```
brew install eigen nlohmann-json
```
### Installation of this project
```
git clone https://github.com/fkgw1228/dnadamage2-modified.git
cd dnadamage2-modified
mkdir build
cd build
cmake ..
make
```

## Usage
Users can view the DNA structure model by runding:
```
./dnadamage2_modified
```
To change which DNA structure is displyed, edit the file path in
`/det/setDNAFile` in `init_vis.mac`. 
To run the application with a macro file, use:
```
./dnadamage2_modified dnadamage2.mac
```
A radom seed can be changed by providing a second argument when running the program.

## Efficiency Tips
The Current geometry is hierarchally defined as follows
- World: the simulation world volume (box)
- Envelope: a volume used for mesuring LET and chemical species yields, and for placing hydrogen atoms in the 'StrandBreakScorer' (orb)
- VoxelStraight: a volume used to contain and replicate the DNA model (box)
At present, the sensitive detector (SD) is attached to **all logical volumes** to score quantities described above. However, this approach can lead to significant computational overhead as the system size increases.
To eveluate DNA damage more efficiently, it is sufficient to attach the SD only to the logical volumes of the DNA model, excluding the volumes listed above. This can be implemented in `DetectorConstruction::ConstructSDAndField()`. Such a modification can substantially improve computational efficiency (speed-up greater than x10!).
