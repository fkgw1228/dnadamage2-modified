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

#include <utility>
#include "G4DNAChemistryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Ellipsoid.hh"
#include "G4DisplacedSolid.hh"
#include "Randomize.hh"
#include "G4UnionSolid.hh"

#include "DNAHandler.hh"
#include "GeometryHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LogicalVolume* DNAGeometryConstructor::CreateDNAGeometry(
    DNAStructure dnaStructure,
    G4bool checkOverlaps
) {
    // Materials
    G4NistManager * man = G4NistManager::Instance();
    G4Material* envelopeWater = man->FindOrBuildMaterial("G4_WATER");
    G4Material* water = man->BuildMaterialWithNewDensity("G4_WATER_COMPOUND",
                                                        "G4_WATER", 1.0*g/cm3);

    DNA dna = dnaStructure.GetDNA();
    // get DNA range
    std::pair<G4ThreeVector, G4ThreeVector> dnaRange = DNAHandler::GetRange(dna);
    G4ThreeVector maxVec = dnaRange.second;
    G4ThreeVector minVec = dnaRange.first;
    // DNA size
    G4double xWidth = (maxVec.x()-minVec.x());
    G4double yWidth = (maxVec.y()-minVec.y());
    G4double zWidth = (maxVec.z()-minVec.z());

    // Envelope solid
    G4String boxNameSolid = "Sol_" + fGeometryName;
    G4Box* boxSolid = new G4Box(boxNameSolid, xWidth/2.0+10*angstrom,
                                              yWidth/2.0+10*angstrom,
                                              zWidth/2.0+10*angstrom);
    // Envelope logic
    G4String boxNameLogic = "Log_" + fGeometryName;
    G4LogicalVolume* boxLogic  = new G4LogicalVolume(boxSolid, 
                                                      envelopeWater, 
                                                      boxNameLogic);

    // Generate a plane
    // lower surface of the box will fit to xy plane
    G4double d = 20 * angstrom;
    G4Box* planeBox = new G4Box("Sol_Plane_Box", d*0.5, d*0.5, d*0.5);
    G4ThreeVector planeDownTranslation(0, 0, -d*0.5);

    // apply translation to the box (xy plane at z=0)
    G4DisplacedSolid* xyPlaneSolid = new G4DisplacedSolid("Sol_XY_Plane_Box",
                                                          planeBox,
                                                          nullptr,
                                                          planeDownTranslation
                                                          );
    size_t ellipsoidCount = 0;

    // loop over all compounds in the DNA structure
    for (auto& strandPair : dna.GetStrands()) {
        G4String chainId = strandPair.first;
        Strand strand = strandPair.second;
        for (auto& nucleotidePair : strand.GetNucleotides()) {
            G4int nucleotideId = nucleotidePair.first;
            Nucleotide nucleotide = nucleotidePair.second;
            for (auto& compoundPair : nucleotide.GetCompounds()) {
                // get compound info
                G4String compoundName = compoundPair.first;
                Compound compound = compoundPair.second;

                // get ellipsoid and planes
                Ellipsoid ellipsoid = dnaStructure.GetEllipsoid(chainId, nucleotideId, compoundName);
                if (ellipsoid.IsEmpty()) continue;  // skip if ellipsoid or planes not found
                // make a solid of the ellipsoid representing the compound
                G4ThreeVector center = ellipsoid.GetCenter();
                G4ThreeVector semiAxisLengths = ellipsoid.GetSemiAxisLengths();
                G4RotationMatrix* axisVectors = new G4RotationMatrix(ellipsoid.GetAxisDirections());

                G4String locationName = chainId + "_" + std::to_string(nucleotideId) + "_" + compoundName;
                // ellipsoid solid of the compound
                G4Ellipsoid* SolidEllipsoid = new G4Ellipsoid("Sol_Ellipsoid_" + locationName,
                                                              semiAxisLengths.x(),
                                                              semiAxisLengths.y(),
                                                              semiAxisLengths.z()
                                                              );
                // apply rotation and translation to the ellipsoid
                G4DisplacedSolid* transformedEllipsoidSolid = 
                                    new G4DisplacedSolid("Sol_Transformed_Ellipsoid_" + locationName,    // Name of the solid
                                                        SolidEllipsoid,                                  // ellipsoid solid      
                                                        axisVectors,                                     // rotation (pointer)
                                                        center                                           // translation
                                                        );
                // ellipsoid to which subtraction will be applied
                G4VSolid* finalSolid = transformedEllipsoidSolid;

                std::vector<Plane> planes = dnaStructure.GetPlanes(chainId, nucleotideId, compoundName);
                size_t planeCount = 0;

                for (Plane& plane : planes) {
                    //if (planeCount != 0) continue;
                    // get rotation matrix and translation vector for the plane
                    G4RotationMatrix* planeRotation = GeometryHandler::GetPlaneRotationMatrix(plane);
                    G4ThreeVector planeTranslation = GeometryHandler::GetPlaneTranslationForEllipsoid(ellipsoid, plane);
                    // apply rotaion and translation to xy plane
                    G4DisplacedSolid* planeSolid = new G4DisplacedSolid("Sol_Transformed_Plane_"+locationName+"_"+std::to_string(planeCount),
                                                                        xyPlaneSolid,
                                                                        planeRotation,
                                                                        planeTranslation
                                                                        );
                    // subtract plane from ellipsoid 
                    finalSolid = new G4SubtractionSolid("Sol_Cut_Ellipsoid_"+locationName+"_"+std::to_string(planeCount++),
                                                        finalSolid,
                                                        planeSolid, nullptr, G4ThreeVector());
                }
                G4LogicalVolume* logicEllipsoid = new G4LogicalVolume(finalSolid,
                                                                      water,
                                                                      "Log_Ellipsoid_" + locationName);
                // place subtracted compound
                G4Colour colour;
                if (G4StrUtil::contains(compoundName, "deoxyribose")) colour = G4Colour::Red();
                else if (G4StrUtil::contains(compoundName, "phosphate")) colour = G4Colour::Yellow();
                else if (G4StrUtil::contains(compoundName, "base")) colour = G4Colour::Green();
                else colour = G4Colour::White();
                G4VisAttributes* MyVisAtt = new G4VisAttributes(colour);
                logicEllipsoid->SetVisAttributes(MyVisAtt);

                G4String ellipsoidName = "Phy_Ellipsoid_" + locationName;

                new G4PVPlacement(0,                // rotation
                                  G4ThreeVector(),  // translation
                                  logicEllipsoid,   // logical volume
                                  ellipsoidName,    // name
                                  boxLogic,        // mother volume
                                  false,            // false
                                  ellipsoidCount++, // copy number
                                  checkOverlaps);           // overlap checking
            }
        }
    }

    return boxLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....