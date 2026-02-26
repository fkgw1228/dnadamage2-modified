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

#include "GeometryHandler.hh"

#include "G4INCLGlobals.hh"
#include "G4SystemOfUnits.hh"

#include "Eigen/Dense"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix* GeometryHandler::GetPlaneRotationMatrix(Plane plane)
{
  /* Obtain the rotation matrix that transforms the xy plane at z=0 to the specified plane */
  G4ThreeVector w = plane.GetW().unit();
  // rotate xy plane and conform upper side of box with true plane
  G4ThreeVector n(0.0, 0.0, 1.0);  // normal vector of xy plane
  G4double dotProduct = n.dot(w);
  G4RotationMatrix* planeRotation = new G4RotationMatrix();
  if (std::abs(dotProduct - 1.0) < 1e-10)
  {
    return planeRotation;  // no rotation needed for parallel planes
  }
  else if (std::abs(dotProduct + 1.0) < 1e-10)
  {
    planeRotation->rotateX(CLHEP::pi);  // 180 degree rotation for anti-parallel planes
    return planeRotation;
  }
  G4double theta = G4INCL::Math::arcCos(dotProduct);  // angle between the plane and xy plane
  // rotation axis
  G4ThreeVector rotationAxis = (n.cross(w)).unit();  // cross product
  G4double n1 = rotationAxis.x(), n2 = rotationAxis.y(), n3 = rotationAxis.z();

  // Rodrigues' rotation formula
  planeRotation->setRows(G4ThreeVector(n1 * n1 * (1 - cos(theta)) + cos(theta),
                                       n1 * n2 * (1 - cos(theta)) + n3 * sin(theta),
                                       n1 * n3 * (1 - cos(theta)) - n2 * sin(theta)),
                         G4ThreeVector(n2 * n1 * (1 - cos(theta)) - n3 * sin(theta),
                                       n2 * n2 * (1 - cos(theta)) + cos(theta),
                                       n2 * n3 * (1 - cos(theta)) + n1 * sin(theta)),
                         G4ThreeVector(n1 * n3 * (1 - cos(theta)) + n2 * sin(theta),
                                       n2 * n3 * (1 - cos(theta)) - n1 * sin(theta),
                                       n3 * n3 * (1 - cos(theta)) + cos(theta)));
  return planeRotation;
}

G4ThreeVector GeometryHandler::GetPlaneTranslationForEllipsoid(Ellipsoid ellipsoid, Plane plane)
{
  /* Return the translation vector from the origin to the center of the intersection
     between the ellipsoid and the plane, The vector is given by:
     c - (b + wc) / (wE^-1w) * (E^-1 w)
     where c is the ellipsoid center, b is the plane bias, w is the plane normal vector,
     and E is the ellipsoid matrix.
  */
  Eigen::Matrix3d E = ellipsoid.GetEMatrix();
  G4ThreeVector c = ellipsoid.GetCenter();
  G4ThreeVector w = plane.GetW().unit();
  G4double b = plane.GetB() / plane.GetW().mag();
  Eigen::Vector3d wEigen;
  wEigen << w.x(), w.y(), w.z();
  Eigen::Vector3d EinW = E.ldlt().solve(wEigen);
  G4double wEin = wEigen.transpose() * EinW;
  G4double scale = (b + w.dot(c)) / wEin;
  G4ThreeVector translation = c - scale * G4ThreeVector(EinW(0), EinW(1), EinW(2));
  return translation;
}