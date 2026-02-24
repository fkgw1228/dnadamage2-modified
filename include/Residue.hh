#ifndef RESIDUE_HH
#define RESIDUE_HH

#include "globals.hh"

#include "Compound.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Residue
{
  public:
    // Constructor
    Residue() = default;
    //Residue(G4int residueId) : fId(residueId), fName("") {}
    Residue(G4int residueId, G4String name) : fId(residueId), fName(name) {}

    // Getter
    G4int GetId() { return this->fId; }
    G4String GetName() { return this->fName; }
    std::map<G4String, Compound>& GetCompounds() { return this->fCompounds; }
    Compound& GetCompound(G4String compoundName);

    // Setter
    void SetName(G4String name) { this->fName = name; }
    void AddCompound(Compound compound) { this->fCompounds[compound.GetName()] = compound; }

  private:
    G4int fId;                                // residue number
    G4String fName;                           // residue name (e.g., DA,DC,DC,DT,)
    std::map<G4String, Compound> fCompounds;  // pairs of name and compound constituting the residue
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif