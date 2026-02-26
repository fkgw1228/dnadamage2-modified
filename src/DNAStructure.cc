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

#include "DNAStructure.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid DNAStructure::GetEllipsoid(G4String chainId, G4int residueId, G4String compoundName)
{
  if (fEllipsoidMap.find(DNAComponentKey{chainId, residueId, compoundName}) == fEllipsoidMap.end())
  {
    return Ellipsoid();  // return empty ellipsoid if not found
  }
  // get ellipsoid from the map
  return fEllipsoidMap.at(DNAComponentKey{chainId, residueId, compoundName});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<Plane> DNAStructure::GetPlanes(G4String chainId, G4int residueId,
                                        G4String compoundName)
{
  if (fPlanesMap.find(DNAComponentKey{chainId, residueId, compoundName}) == fPlanesMap.end())
  {
    return std::vector<Plane>();  // return empty vector if not found
  }
  // get planes from the map
  return fPlanesMap.at(DNAComponentKey{chainId, residueId, compoundName});
}

void DNAStructure::SetEllipsoid(G4String chainId, G4int residueId, G4String compoundName,
                             Ellipsoid ellipsoid)
{
  fEllipsoidMap[DNAComponentKey{chainId, residueId, compoundName}] = ellipsoid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::AddPlane(G4String chainId, G4int residueId, G4String compoundName, Plane plane)
{
  fPlanesMap[DNAComponentKey{chainId, residueId, compoundName}].push_back(plane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::SetPlanes(G4String chainId, G4int residueId, G4String compoundName,
                          std::vector<Plane> planes)
{
  fPlanesMap[DNAComponentKey{chainId, residueId, compoundName}] = planes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......