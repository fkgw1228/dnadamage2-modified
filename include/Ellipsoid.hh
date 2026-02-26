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
//
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

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