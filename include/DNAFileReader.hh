#ifndef DNAREADER_HH
#define DNAREADER_HH

#include "G4ThreeVector.hh"

#include "Atom.hh"
#include "Chain.hh"
#include "Compound.hh"
#include "Ellipsoid.hh"
#include "Model.hh"
#include "Plane.hh"
#include "Residue.hh"
#include "DNAStructure.hh"
#include <nlohmann/json.hpp>

#include <string>
#include <vector>

using json = nlohmann::json;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAFileReader
{
  public:
    // Constructor
    DNAFileReader(G4String fileName) : fFileName(fileName) {}

    // Setter
    void SetFileName(G4String fileName) { this->fFileName = fileName; }

    DNAStructure ReadDNAFile(G4int modelId);

  private:
    Model GenerateModel(G4int modelId, const json& dnaInfo);
    Chain GenerateChain(G4String chainId, const json& chainInfo);
    Residue GenerateResidue(G4int residueId, const json& residueInfo);
    Compound GenerateCompound(G4String compoundName, const json& compoundInfo);
    Atom GenerateAtom(const json& atomInfo);
    DNAStructure GenerateStructure(const json& dnaInfo);
    Ellipsoid GenerateEllipsoid(const json& ellipsoidInfo);
    Plane GeneratePlane(const json& planeInfo);

    G4String fFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif