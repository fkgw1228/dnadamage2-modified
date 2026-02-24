#ifndef CHAIN_HH
#define CHAIN_HH

#include "globals.hh"

#include "Residue.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Chain
{
  public:
    // Constructor
    Chain() = default;
    Chain(G4String id) : fId(id) {}

    // Getter
    G4String GetId() { return fId; }
    std::map<G4int, Residue>& GetResidues() { return fResidues; }
    Residue& GetResidue(G4int residueId);

    // Setter
    void AddResidue(Residue residue) { fResidues[residue.GetId()] = residue; }

  private:
    G4String fId;                       // Chain id
    std::map<G4int, Residue> fResidues; // Constituent residues
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif