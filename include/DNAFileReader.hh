#ifndef DNAREADER_HH
#define DNAREADER_HH

#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "G4ThreeVector.hh"
#include "DNAStructure.hh"
#include "DNA.hh"
#include "Strand.hh"
#include "Nucleotide.hh"
#include "Compound.hh"
#include "Atom.hh"
#include "Ellipsoid.hh"
#include "Plane.hh"

using json = nlohmann::json;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAFileReader {
public:
    // Constructor
    DNAFileReader(G4String fileName) : fFileName(fileName) {}
    
    // Setter
    void SetFileName(G4String fileName) {this->fFileName = fileName;}

    DNAStructure ReadDNAFile(G4int dnaId);

private:
    DNA          GenerateDNA(G4int dnaId, const json& dnaInfo);
    Strand       GenerateStrand(G4String chainId, const json& strandInfo);
    Nucleotide   GenerateNucleotide(G4int nucleotideId, const json& nucleotideInfo);
    Compound     GenerateCompound(G4String compoundName, const json& compoundInfo);
    Atom         GenerateAtom(const json& atomInfo);
    DNAStructure GenerateGeometry(const json& dnaInfo);
    Ellipsoid    GenerateEllipsoid(const json& ellipsoidInfo);
    Plane        GeneratePlane(const json& planeInfo);

    G4String fFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif