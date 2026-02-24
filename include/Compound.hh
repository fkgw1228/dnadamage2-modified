#ifndef COMPOUND_HH
#define COMPOUND_HH

#include "globals.hh"

#include "Atom.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Compound
{
  public:
    // Constructor
    Compound() = default;
    Compound(G4String name) : fName(name) {}

    // Getters
    G4String GetName() { return fName; }
    std::vector<Atom>& GetAtoms() { return fAtoms; }

    // Setters
    void AddAtom(Atom atom) { this->fAtoms.push_back(atom); }
    void SetAtoms(std::vector<Atom> atoms) { this->fAtoms = atoms; }

  private:
    G4String fName;            // compound name
    std::vector<Atom> fAtoms;  // constituent atom list
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif