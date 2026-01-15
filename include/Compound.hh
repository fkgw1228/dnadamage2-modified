#ifndef COMPOUND_HH
#define COMPOUND_HH

#include <vector>
#include "globals.hh"
#include "Atom.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Compound {
public:
    // Constructor
    Compound() : fName("") {}
    Compound(G4String name): fName(name) {}

    // Getters
    G4String GetName() const { return fName; }
    const std::vector<Atom>& GetAtoms() const { return fAtoms; }

    // Setters
    void AddAtom(Atom atom) { this->fAtoms.push_back(atom); }
    void SetAtoms(std::vector<Atom> atoms) { this->fAtoms = atoms; }

private:
    G4String fName;                 // compound name
    std::vector<Atom> fAtoms;       // atom list
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif