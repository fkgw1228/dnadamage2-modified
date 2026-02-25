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
// Authors: J. Naoki D. Kondo (UCSF, US) : 10/10/2021
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file DNAGeometryConstructor.cc
/// \brief Implementation of the plasmid load methods for the geometry

#include "DNAGeometryConstructor.hh"

#include "G4Box.hh"
#include "G4DisplacedSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"

#include "DNAHandler.hh"
#include "GeometryHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LogicalVolume* DNAGeometryConstructor::CreateDNAGeometry(DNAStructure& structure)
{
  // Materials
  // * Material for the compound volumes has different name for chemistry control (see
  // TimeStepAction)
  G4NistManager* man = G4NistManager::Instance();
  G4Material* envelopeWater = man->FindOrBuildMaterial("G4_WATER");
  G4Material* water =
    man->BuildMaterialWithNewDensity("G4_WATER_COMPOUND", "G4_WATER", 1.0 * g / cm3);

  Model model = structure.GetModel();
  auto [minEdge, maxEdge] = DNAHandler::GetRange(model);
  G4ThreeVector voxelSize = maxEdge - minEdge
                              + G4ThreeVector(50*angstrom,50*angstrom,50*angstrom);

  G4LogicalVolume* boxLogic = CreateVoxelBox(voxelSize, envelopeWater);
  G4DisplacedSolid* xyPlaneSolid = BuildXYPlane(20*angstrom);

  size_t ellipsoidCount = 0;

  // loop over all compounds in the DNA structure
  for (auto& [chainId, chain] : model.GetChains())
  {
    for (auto& [residueId, residue] : chain.GetResidues())
    {
      for (auto& [compoundName, compound] : residue.GetCompounds())
      {
        // get ellipsoid and planes
        Ellipsoid ellipsoid = structure.GetEllipsoid(chainId, residueId, compoundName);
        std::vector<Plane> planes = structure.GetPlanes(chainId, residueId, compoundName);
        if (ellipsoid.IsEmpty()) continue;  // skip if ellipsoid not found

        G4VSolid* compoundSolid = 
                BuildCompoundSolid(
                  chainId, residueId, compoundName,
                  ellipsoid, planes, xyPlaneSolid
                );
        
        PlaceComound(boxLogic, compoundSolid,
                    chainId, residueId, compoundName,
                    water, ellipsoidCount++);
      }
    }
  }

  return boxLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LogicalVolume* DNAGeometryConstructor::CreateVoxelBox(G4ThreeVector voxelSize,
                                                        G4Material* envelopeWater)
{
  // Envelope solid
  G4String boxNameSolid = "Sol_" + fGeometryName;
  G4Box* boxSolid = new G4Box(boxNameSolid, voxelSize.x()*0.5, voxelSize.y()*0.5, voxelSize.z()*0.5);
  // Envelope logic
  G4String boxNameLogic = "Log_" + fGeometryName;
  G4LogicalVolume* boxLogic = new G4LogicalVolume(boxSolid, envelopeWater, boxNameLogic);

  return boxLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4DisplacedSolid* DNAGeometryConstructor::BuildXYPlane(G4double size)
{
  // upper surface of the box will fit to xy plane
  G4Box* planeBox = new G4Box("Sol_Plane_Box", size * 0.5, size * 0.5, size * 0.5);
  G4ThreeVector planeDownTranslation(0, 0, -size * 0.5);
  // apply translation to the box (xy plane at z=0)
  return new G4DisplacedSolid("Sol_XY_Plane_Box", planeBox, nullptr, planeDownTranslation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VSolid* DNAGeometryConstructor::BuildCompoundSolid(
            G4String chainId, G4int residueId, G4String compoundName,
            Ellipsoid ellipsoid, std::vector<Plane> planes,
            G4DisplacedSolid* xyPlaneSolid)
{
  G4String locationName = chainId + "_" + std::to_string(residueId) + "_" + compoundName;

  // make a solid of the ellipsoid representing the compound
  G4ThreeVector center = ellipsoid.GetCenter();
  G4ThreeVector semiAxisLengths = ellipsoid.GetSemiAxisLengths();
  G4RotationMatrix* axisDirections = new G4RotationMatrix(ellipsoid.GetAxisDirections());
  
  G4Ellipsoid* SolidEllipsoid =
      new G4Ellipsoid("Sol_Ellipsoid_" + locationName,
                      semiAxisLengths.x(),
                      semiAxisLengths.y(),
                      semiAxisLengths.z());
  
  // First, apply only rotation to the ellipsoid (centered at origin)
  G4DisplacedSolid* transformedEllipsoidSolid =
    new G4DisplacedSolid("Sol_Transformed_Ellipsoid_" + locationName,
                         SolidEllipsoid,    // ellipsoid solid
                         axisDirections,    // rotation
                         center);           // translation

  G4VSolid* cutEllipsoidSolid = ApplyPlaneCuts(
      transformedEllipsoidSolid,
      ellipsoid, planes, locationName, xyPlaneSolid
    );

  return cutEllipsoidSolid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VSolid* DNAGeometryConstructor::ApplyPlaneCuts(
              G4DisplacedSolid* ellipsoidSolid,
              Ellipsoid ellipsoid,
              std::vector<Plane> planes,
              G4String locationName,
              G4DisplacedSolid* xyPlaneSolid)
{
  size_t planeCount = 0;
  // ellipsoid to which subtraction will be applied (still at origin)
  G4VSolid* finalSolid = ellipsoidSolid;

  for (Plane& plane : planes)
  {
    // get rotation matrix and translation vector for the plane
    G4RotationMatrix* planeRotation = GeometryHandler::GetPlaneRotationMatrix(plane);
    G4ThreeVector planeTranslation =
      GeometryHandler::GetPlaneTranslationForEllipsoid(ellipsoid, plane);

    // apply rotation and translation to xy plane
    G4DisplacedSolid* planeSolid = new G4DisplacedSolid(
      "Sol_Transformed_Plane_" + locationName + "_" + std::to_string(planeCount),
      xyPlaneSolid, planeRotation, planeTranslation);
    // subtract plane from ellipsoid
    finalSolid = new G4SubtractionSolid("Sol_Cut_Ellipsoid_" + locationName + "_"
                                          + std::to_string(planeCount++),
                                        finalSolid, planeSolid, nullptr, G4ThreeVector());
  }

  return finalSolid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DNAGeometryConstructor::PlaceComound(
              G4LogicalVolume* mother,
              G4VSolid* solid,
              G4String chainId, G4int residueId, G4String compoundName,
              G4Material* material,
              G4int copyNumber
            )
{
  G4String locationName = chainId + "_" + std::to_string(residueId) + "_" + compoundName;
  G4LogicalVolume* logicEllipsoid =
          new G4LogicalVolume(solid, material, "Log_Ellipsoid_" + locationName);
  // place subtracted compound
  G4VisAttributes* visAttr = fVisAttrMap[compoundName];
  logicEllipsoid->SetVisAttributes(visAttr);

  G4String ellipsoidName = "Phy_Ellipsoid_" + locationName;

  // Place the ellipsoid
  new G4PVPlacement(0,                // no rotation
                    G4ThreeVector(),  // no translation
                    logicEllipsoid,   // logical volume
                    ellipsoidName,    // name
                    mother,           // mother volume
                    false,            // false
                    copyNumber,       // copy number
                    fCheckOverlaps);  // overlap checking
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// G4LogicalVolume* DNAGeometryConstructor::CreateDNAGeometry(DNAStructure& structure,
//                                                            G4double voxelSize)
// {
//   // Materials
//   // * Material for the compound volumes has different name for chemistry control (see
//   // TimeStepAction)
//   G4NistManager* man = G4NistManager::Instance();
//   G4Material* envelopeWater = man->FindOrBuildMaterial("G4_WATER");
//   G4Material* water =
//     man->BuildMaterialWithNewDensity("G4_WATER_COMPOUND", "G4_WATER", 1.0 * g / cm3);

//   Model model = structure.GetModel();

//   // Envelope solid
//   G4String boxNameSolid = "Sol_" + fGeometryName;
//   G4Box* boxSolid = new G4Box(boxNameSolid, voxelSize * 0.5, voxelSize * 0.5, voxelSize * 0.5);
//   // Envelope logic
//   G4String boxNameLogic = "Log_" + fGeometryName;
//   G4LogicalVolume* boxLogic = new G4LogicalVolume(boxSolid, envelopeWater, boxNameLogic);

//   // Generate a plane
//   // upper surface of the box will fit to xy plane
//   G4double d = 20 * angstrom;
//   G4Box* planeBox = new G4Box("Sol_Plane_Box", d * 0.5, d * 0.5, d * 0.5);
//   G4ThreeVector planeDownTranslation(0, 0, -d * 0.5);
//   // apply translation to the box (xy plane at z=0)
//   G4DisplacedSolid* xyPlaneSolid =
//     new G4DisplacedSolid("Sol_XY_Plane_Box", planeBox, nullptr, planeDownTranslation);
//   size_t ellipsoidCount = 0;

//   // loop over all compounds in the DNA structure
//   for (auto& [chainId, chain] : model.GetChains())
//   {
//     for (auto& [residueId, residue] : chain.GetResidues())
//     {
//       for (auto& [compoundName, compound] : residue.GetCompounds())
//       {
//         // get ellipsoid and planes
//         Ellipsoid ellipsoid = structure.GetEllipsoid(chainId, residueId, compoundName);
//         if (ellipsoid.IsEmpty()) continue;  // skip if ellipsoid or planes not found
//         // make a solid of the ellipsoid representing the compound
//         G4ThreeVector center = ellipsoid.GetCenter();
//         G4ThreeVector semiAxisLengths = ellipsoid.GetSemiAxisLengths();
//         G4RotationMatrix axisDirections = ellipsoid.GetAxisDirections();

//         G4String locationName = chainId + "_" + std::to_string(residueId) + "_" + compoundName;
//         // ellipsoid solid of the compound (at origin, aligned with axes)
//         G4Ellipsoid* SolidEllipsoid =
//           new G4Ellipsoid("Sol_Ellipsoid_" + locationName, semiAxisLengths.x(), semiAxisLengths.y(),
//                           semiAxisLengths.z());
//         // First, apply only rotation to the ellipsoid (centered at origin)
//         G4RotationMatrix* axisVectors = new G4RotationMatrix(axisDirections);
//         G4DisplacedSolid* transformedEllipsoidSolid =
//           new G4DisplacedSolid("Sol_Transformed_Ellipsoid_" + locationName,  // Name of the solid
//                                SolidEllipsoid,  // ellipsoid solid
//                                axisVectors,  // rotation (pointer)
//                                center  // translation
//           );
//         // Planes cutting the ellipsoid
//         std::vector<Plane> planes = structure.GetPlanes(chainId, residueId, compoundName);
//         size_t planeCount = 0;
//         // ellipsoid to which subtraction will be applied (still at origin)
//         G4VSolid* finalSolid = transformedEllipsoidSolid;

//         for (Plane& plane : planes)
//         {
//           // get rotation matrix and translation vector for the plane
//           G4RotationMatrix* planeRotation = GeometryHandler::GetPlaneRotationMatrix(plane);
//           G4ThreeVector planeTranslation =
//             GeometryHandler::GetPlaneTranslationForEllipsoid(ellipsoid, plane);

//           // apply rotation and translation to xy plane
//           G4DisplacedSolid* planeSolid = new G4DisplacedSolid(
//             "Sol_Transformed_Plane_" + locationName + "_" + std::to_string(planeCount),
//             xyPlaneSolid, planeRotation, planeTranslation);
//           // subtract plane from ellipsoid
//           finalSolid = new G4SubtractionSolid("Sol_Cut_Ellipsoid_" + locationName + "_"
//                                                 + std::to_string(planeCount++),
//                                               finalSolid, planeSolid, nullptr, G4ThreeVector());
//         }
//         G4LogicalVolume* logicEllipsoid =
//           new G4LogicalVolume(finalSolid, water, "Log_Ellipsoid_" + locationName);
//         // place subtracted compound
//         G4VisAttributes* visAttr = fVisAttrMap[compoundName];
//         logicEllipsoid->SetVisAttributes(visAttr);

//         G4String ellipsoidName = "Phy_Ellipsoid_" + locationName;

//         // Place the ellipsoid
//         new G4PVPlacement(0,  // no rotation
//                           G4ThreeVector(),  // no translation
//                           logicEllipsoid,  // logical volume
//                           ellipsoidName,  // name
//                           boxLogic,  // mother volume
//                           false,  // false
//                           ellipsoidCount++,  // copy number
//                           fCheckOverlaps);  // overlap checking
//       }
//     }
//   }

//   return boxLogic;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
