#ifndef PLANE_HH
#define PLANE_HH

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Plane {
public:
    // Constructor
    Plane(G4ThreeVector w, G4double b): w(w), b(b) {}

    // Getter
    G4ThreeVector GetW() const { return w; }
    G4double GetB() const { return b; }

    // Setter
    void SetW(G4ThreeVector w) { this->w = w; }
    void SetB(G4double b) { this->b = b; }

private:
    G4ThreeVector w;  // normal vector
    G4double b;       // bias
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif