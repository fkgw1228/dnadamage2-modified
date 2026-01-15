#ifndef DNASTRUCTURE_HH
#define DNASTRUCTURE_HH

#include <vector>
#include <map>
#include "DNA.hh"
#include "Ellipsoid.hh"
#include "Plane.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct DNALocationKey {
    G4String chainId;
    G4int nucleotideId;
    G4String compoundName;

    // Define less-than operator for use in std::map
    bool operator<(const DNALocationKey& other) const {
        if (chainId != other.chainId)             return chainId < other.chainId;
        if (nucleotideId != other.nucleotideId)   return nucleotideId < other.nucleotideId;
        return compoundName < other.compoundName;
    }
};

typedef std::map<DNALocationKey, Ellipsoid> DNAEllipsoidMap;
typedef std::map<DNALocationKey, std::vector<Plane>> DNAPlanesMap;

class DNAStructure {
public:
    // Constructor
    DNAStructure() {}
    DNAStructure(DNA dna) : fDna(dna) {}

    // Getter
    const DNA& GetDNA() const { return fDna; }
    Ellipsoid GetEllipsoid(G4String chainId, G4int nucleotideId, G4String compoundName) const;
    std::vector<Plane> GetPlanes(G4String chainId, G4int nucleotideId, G4String compoundName) const;

    // Setter
    void SetDNA(DNA dna) { fDna = dna; }
    void SetEllipsoid(G4String chainId, G4int nucleotideId, G4String compoundName, Ellipsoid ellipsoid);
    void AddPlane(G4String chainId, G4int nucleotideId, G4String compoundName, Plane plane);
    void SetPlanes(G4String chainId, G4int nucleotideId, G4String compoundName, std::vector<Plane> planes);

private:
    DNA fDna;
    DNAEllipsoidMap fEllipsoidMap;
    DNAPlanesMap fPlanesMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

# endif