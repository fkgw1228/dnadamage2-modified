#ifndef ELLIPSOID_HH
#define ELLIPSOID_HH

#include <iostream>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include <Eigen/Dense>

using Matrix3d = Eigen::Matrix3d;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Ellipsoid {
public:
    // Constructor
    Ellipsoid ()
        : fEMatrix(Matrix3d::Zero()),
            fCenter(G4ThreeVector(0,0,0)),
            fSemiAxisLengths(G4ThreeVector(0,0,0)),
            fAxisDirections(G4RotationMatrix()),
            fIsEmpty(true) {}
    Ellipsoid(Matrix3d E_matrix, G4ThreeVector center);

    // Getter
    Matrix3d GetEMatrix() const { return fEMatrix; }
    G4ThreeVector GetCenter() const { return fCenter; }
    G4ThreeVector GetSemiAxisLengths() const { return fSemiAxisLengths; }
    G4RotationMatrix GetAxisDirections() const { return fAxisDirections; }
    
    // setter
    void SetCenter(G4ThreeVector center) {this->fCenter = center;}

    G4bool IsEmpty() const { return fIsEmpty; };

private:
    Matrix3d fEMatrix;                  // matrix representing the shape of ellipsoid (3,3)
    G4ThreeVector fCenter;              // ellipsoid center
    G4ThreeVector fSemiAxisLengths;     // ellipsoid axis vector (3,)
    G4RotationMatrix fAxisDirections;   // ellipsoid axis directions (3,3)
    G4bool fIsEmpty;                    // flag to indicate if the ellipsoid is empty
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif