//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

#include "Ellipsoid.hh"

#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector Ellipsoid::GetSemiAxisLengths()
{
  // calculate semi-axis lengths and axis directions from E matrix
  Eigen::SelfAdjointEigenSolver<Matrix3d> solver(fEMatrix);
  Eigen::Vector3d eigenvalues = solver.eigenvalues();

  // semi-axis lengths are given by 1/sqrt(eigenvalue)
  return  G4ThreeVector(1.0 / sqrt(eigenvalues(0)) * angstrom,
                        1.0 / sqrt(eigenvalues(1)) * angstrom,
                        1.0 / sqrt(eigenvalues(2)) * angstrom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix Ellipsoid::GetAxisDirections()
{
  // Calculate axis directions from E matrix
  Eigen::SelfAdjointEigenSolver<Matrix3d> solver(fEMatrix);

  // axis directions are given by eigenvectors
  Matrix3d eigenvecs = solver.eigenvectors();
  eigenvecs.col(0).normalize();
  eigenvecs.col(2) = (eigenvecs.col(0).cross(eigenvecs.col(1))).normalized();
  eigenvecs.col(1) = (eigenvecs.col(2).cross(eigenvecs.col(0))).normalized();

  // G4AffineTransform applies vectors as v*R, therefore R must be V^T.
  G4RotationMatrix axisDirections;
  axisDirections.setRows(
    G4ThreeVector(eigenvecs(0, 0), eigenvecs(1, 0), eigenvecs(2, 0)),
    G4ThreeVector(eigenvecs(0, 1), eigenvecs(1, 1), eigenvecs(2, 1)),
    G4ThreeVector(eigenvecs(0, 2), eigenvecs(1, 2), eigenvecs(2, 2))
  );
  return axisDirections;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
