#ifndef ATOM_HH
#define ATOM_HH

#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Atom
{
  public:
    // Constructor
    Atom(G4String name, G4int id, G4String element, G4ThreeVector position, G4double radius)
      : fName(name), fId(id), fElement(element), fPosition(position), fRadius(radius) {}

    // Getters
    G4String GetName() { return fName; }
    G4int GetId() { return fId; }
    G4String GetElement() { return fElement; }
    G4ThreeVector GetPosition() { return fPosition; }
    G4double GetRadius() { return fRadius; }

    // Setter
    void SetPosition(G4ThreeVector newPosition) { fPosition = newPosition; }

  private:
    G4String fName;           // Atom name
    G4int fId;                // Atom ID
    G4String fElement;        // Element
    G4ThreeVector fPosition;  // Atomic position
    G4double fRadius;         // Atom's radius
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif