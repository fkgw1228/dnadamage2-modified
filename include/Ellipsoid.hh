#ifndef ELLIPSOID_HH
#define ELLIPSOID_HH

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include <Eigen/Dense>

#include <iostream>

using Matrix3d = Eigen::Matrix3d;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Ellipsoid
{
  public:
    //constructor
    Ellipsoid(): fIsEmpty(true) {}
    Ellipsoid(Matrix3d EMatrix, G4ThreeVector center):
      fEMatrix(EMatrix), fCenter(center), fIsEmpty(false) {}

    // Getter
    Matrix3d GetEMatrix() { return fEMatrix; }
    G4ThreeVector GetCenter() { return fCenter; }
    G4ThreeVector GetSemiAxisLengths();
    G4RotationMatrix GetAxisDirections();

    // setter
    void SetCenter(G4ThreeVector center) { this->fCenter = center; }

    G4bool IsEmpty() { return fIsEmpty; };

  private:
    Matrix3d fEMatrix;                // matrix representing the shape of ellipsoid (3,3)
    G4ThreeVector fCenter;            // ellipsoid center
    G4bool fIsEmpty;                  // flag to indicate if the ellipsoid is empty
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif