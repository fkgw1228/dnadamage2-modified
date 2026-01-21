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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Update: H. Tran (IRSN, France) :20/12/2018
//         J. Naoki D. Kondo (UCSF, US): 10/10/2021
//
/// \file PhysGeoImport.hh
/// \brief Definition of the plasmid load methods for the geometry

#ifndef DNADAMAGE2_GeoImport_h
#define DNADAMAGE2_GeoImport_h 1

#include <fstream>
#include <algorithm>

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include <memory>
#include "G4H2O.hh"
#include "G4Electron_aq.hh"

#include "DNAStructure.hh"

class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAGeometryConstructor
{
public:
    // Constructor
    DNAGeometryConstructor() {};
    ~DNAGeometryConstructor() = default;

    G4LogicalVolume* CreateDNAGeometry(DNAStructure dnaStructure, G4bool checkOverlaps);

private:
    G4String fGeometryName = "Voxel_Straight";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
