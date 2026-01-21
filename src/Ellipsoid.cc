#include "Ellipsoid.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid::Ellipsoid(Matrix3d E_matrix, G4ThreeVector center) 
    : fEMatrix(E_matrix), fCenter(center) {
    // calculate semi-axis lengths and axis directions from E matrix
    Eigen::SelfAdjointEigenSolver<Matrix3d> solver(E_matrix);
    Eigen::Vector3d eigenvalues = solver.eigenvalues();

    // semi-axis lengths are given by 1/sqrt(eigenvalue)
    fSemiAxisLengths = G4ThreeVector(
        1.0 / sqrt(eigenvalues(0)) * angstrom,
        1.0 / sqrt(eigenvalues(1)) * angstrom,
        1.0 / sqrt(eigenvalues(2)) * angstrom
    );

    // axis directions are given by eigenvectors
    Matrix3d eigenvecs = solver.eigenvectors();
    eigenvecs.col(0).normalize();
    eigenvecs.col(2) = (eigenvecs.col(0).cross(eigenvecs.col(1))).normalized();
    eigenvecs.col(1) = (eigenvecs.col(2).cross(eigenvecs.col(0))).normalized();

    // set rotation matrix from eigenvectors
    double mxx = eigenvecs(0,0), mxy = eigenvecs(1,0), mxz = eigenvecs(2,0);
    double myx = eigenvecs(0,1), myy = eigenvecs(1,1), myz = eigenvecs(2,1);
    double mzx = eigenvecs(0,2), mzy = eigenvecs(1,2), mzz = eigenvecs(2,2);

    G4RotationMatrix rot;
    rot.setRows(
        G4ThreeVector(mxx, mxy, mxz),
        G4ThreeVector(myx, myy, myz),
        G4ThreeVector(mzx, mzy, mzz)
    );

    fAxisDirections = rot;
    fIsEmpty = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......