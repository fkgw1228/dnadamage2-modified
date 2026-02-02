README: [Japanese](README.ja.md)

# dnadamage2-modified

## Requirements

## Installation
```
git clone https://github.com/fkgw1228/dnadamage2-modified.git
cd dnadamage2-modified
mkdir build
cd build
cmake ..
make
```
Linux:
```
sudo apt update
sudo apt install -y libeigne3-dev nlohmann-json3-dev
```
Linux by not using apt
```
# A directory for installation
mkdir ~/local
# Eigen3 installation
git clone https://gitlab.com/libeigen/eigen.git
cd eigen && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIFX=~/local
make install
# nlohmann_json installtion
git clone https;//github.com/nlohmann/json.git
cd json && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/local -DJSON_BuildTests=OFF
make install
# Installation check
ls ~/local/include         # eigen3/  nlohmann/
```
Mac:
```
brew install eigen nlohmann-json
```
