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

#ifndef DNAGEOMETRYCONSTRUCTOR_HH
#define DNAGEOMETRYCONSTRUCTOR_HH 1

#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include "Ellipsoid.hh"
#include "Plane.hh"
#include "DNAStructure.hh"

class G4LogicalVolume;
class G4VSolid;
class G4DisplacedSolid;
class G4SubtractionSolid;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAGeometryConstructor
{
  public:
    // Constructor
    DNAGeometryConstructor() = default;
    DNAGeometryConstructor(G4bool checkOverlaps)
        : DNAGeometryConstructor() { fCheckOverlaps = checkOverlaps; }
    ~DNAGeometryConstructor() = default;

    G4LogicalVolume* CreateDNAGeometry(DNAStructure&);

  
  private:
    G4LogicalVolume* CreateVoxelBox(G4ThreeVector voxelSize, G4Material* envelopeWater);
    G4DisplacedSolid* BuildXYPlane(G4double size);
    G4VSolid* BuildCompoundSolid(
              G4String chainId, G4int residueId, G4String compoundName,
              Ellipsoid ellipsoid, std::vector<Plane> planes,
              G4DisplacedSolid* xyPlaneSolid);
    G4VSolid* ApplyPlaneCuts(
              G4DisplacedSolid* ellipsoidSolid,
              Ellipsoid ellipsoid,
              std::vector<Plane> planes,
              G4String locationName,
              G4DisplacedSolid* xyPlaneSolid);
    void PlaceComound(
              G4LogicalVolume* mother,
              G4VSolid* solid,
              G4String chainId, G4int residueId, G4String compoundName,
              G4Material* compoundWater,
              G4int copyNumber
            );

    G4String fGeometryName = "Voxel_Straight";
    G4bool fCheckOverlaps = false;
    std::map<G4String, G4VisAttributes*> fVisAttrMap
            = {
              {"deoxyribose", new G4VisAttributes(G4Color::Red())},
              {"phosphate",   new G4VisAttributes(G4Color::Yellow())},
              {"base",        new G4VisAttributes(G4Color::Green())},
              {"default",     new G4VisAttributes(G4Color::White())}
              };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
