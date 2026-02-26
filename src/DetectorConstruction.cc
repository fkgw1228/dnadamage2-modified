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
// This code is derived from the dnadamage2 example.
// dnadamage2 example authors: J. Naoki D. Kondo (UCSF, US),
//                             J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UnionSolid.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"

#include "DNAGeometryConstructor.hh"
#include "DNAFileReader.hh"
#include "DNAHandler.hh"
#include "ScoreLET.hh"
#include "ScoreSpecies.hh"
#include "ScoreStrandBreaks.hh"

#include <chrono>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {
  fpDetDir = new G4UIdirectory("/det/");
  fpDetDir->SetGuidance("Detector control.");

  fpWorldSizeUI = new G4UIcmdWithADoubleAndUnit("/det/setWorldSize", this);
  fpWorldSizeUI->SetGuidance("Set the full size of the world box.");
  fpWorldSizeUI->SetDefaultUnit("um");

  fpEnvelopeDiameterUI =
      new G4UIcmdWithADoubleAndUnit("/det/setEnvelopeDiameter", this);
  fpEnvelopeDiameterUI->SetGuidance("Set the diameter of the DNA envelope.");
  fpEnvelopeDiameterUI->SetDefaultUnit("um");

  fpOffsetFileUI = new G4UIcmdWithAString("/det/setOffsetFile", this);
  fpOffsetFileUI->SetGuidance("Set the DNA offset file name.");

  fpNumberOfDNAUI = new G4UIcmdWithAnInteger("/det/setNumberOfDNA", this);
  fpNumberOfDNAUI->SetGuidance("Set the number of DNA structures to be placed.");

  fpDNAFileUI = new G4UIcmdWithAString("/det/setDNAFile", this);
  fpDNAFileUI->SetGuidance("Set the DNA structure file name (json).");

  fpEnableDNAVolumesUI = new G4UIcmdWithABool("/det/enableDNAVolumes", this);
  fpEnableDNAVolumesUI->SetGuidance("Enable or disable DNA volumes for simulation.");

  fpCheckOverlapsUI = new G4UIcmdWithABool("/det/checkOverlaps", this);
  fpCheckOverlapsUI->SetGuidance("Enable or disable overlap checking during geometry construction.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction() {
  delete fpDetDir;
  delete fpWorldSizeUI;
  delete fpEnvelopeDiameterUI;
  delete fpOffsetFileUI;
  delete fpDNAFileUI;
  delete fpEnableDNAVolumesUI;
  delete fpNumberOfDNAUI;
  delete fpCheckOverlapsUI;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetNewValue(G4UIcommand *command,
                                       G4String newValue) {
  if (command == fpWorldSizeUI) {
    fWorldSize = fpWorldSizeUI->GetNewDoubleValue(newValue);
  }
  if (command == fpEnvelopeDiameterUI) {
    fEnvelopeDiameter = fpEnvelopeDiameterUI->GetNewDoubleValue(newValue);
  }
  if (command == fpNumberOfDNAUI) {
    fNumberOfDNA = fpNumberOfDNAUI->GetNewIntValue(newValue);
  }
  if (command == fpOffsetFileUI) {
    ReadOffsetFile(G4String(newValue));
  }
  if (command == fpDNAFileUI) {
    fDNAFile = newValue;
  }
  if (command == fpEnableDNAVolumesUI) {
    fEnableDNAVolumes = fpEnableDNAVolumesUI->GetNewBoolValue(newValue);
  }
  if (command == fpCheckOverlapsUI) {
    fCheckOverlaps = fpCheckOverlapsUI->GetNewBoolValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume *DetectorConstruction::Construct() {
  // Water is defined from NIST material database
  // The name of the world material is different for chemistry control (see TimeStepAction)
  G4NistManager *man = G4NistManager::Instance();
  G4Material *normalWater = man->FindOrBuildMaterial("G4_WATER");
  G4Material *worldWater = man->BuildMaterialWithNewDensity(
      "G4_WATER_WORLD", "G4_WATER", 1.0 * g / cm / cm / cm);

  // World
  G4Box *worldBox = new G4Box("Sol_World",
                              0.5 * fWorldSize,
                              0.5 * fWorldSize,
                              0.5 * fWorldSize);
  G4LogicalVolume *logicWorld = 
                    new G4LogicalVolume(worldBox,
                                        worldWater,
                                        "Log_World");
  G4VPhysicalVolume *physWorld =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        "Phy_World",     // name of physical volume
                        logicWorld,      // logical volume
                        0,               // no mother volume
                        false,           // no use
                        0,               // copy number
                        false);          // checking overlaps
  // World Vis Attributes
  G4VisAttributes *worldVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
  worldVis->SetVisibility(false);
  logicWorld->SetVisAttributes(worldVis);

  // Envelope for DNA geometry (sphere)
  G4Orb *dnaEnvelope = new G4Orb("Sol_DNA_Envelope", fEnvelopeDiameter * 0.5);
  G4LogicalVolume *logicDNAEnvelope =
      new G4LogicalVolume(dnaEnvelope, normalWater, "Log_DNA_Envelope");
  new G4PVPlacement(0,                  // no rotation
                    G4ThreeVector(),    // at (0,0,0)
                    logicDNAEnvelope,   // its logical volume
                    "Phy_DNA_Envelope", // its name
                    logicWorld,         // its mother  volume
                    false,              // no use
                    0,                  // copy number
                    false);             // checking overlaps

  // DNA Envelope Vis Attributes
  G4VisAttributes *dnaEnvelopeVis =
      new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
  dnaEnvelopeVis->SetForceSolid(true);
  logicDNAEnvelope->SetVisAttributes(dnaEnvelopeVis);

  // DNA geometry construction
  auto start = std::chrono::high_resolution_clock::now();
  DNAFileReader reader = DNAFileReader(fDNAFile);
  DNAStructure structure = reader.ReadDNAFile(1);
  
  // DNAHandler::Print(structure);

  // Voxel for storing DNA geometry
  DNAGeometryConstructor constructor = DNAGeometryConstructor(fCheckOverlaps);
  G4LogicalVolume *logicStraightVoxel =
      constructor.CreateDNAGeometry(structure);

  for (G4int i = 1; i <= fNumberOfDNA; i++) {
    if (fEnableDNAVolumes)
      new G4PVPlacement(0,                    // no rotation
                        fOffsets[i - 1],      // at offset position
                        logicStraightVoxel,   // its logical volume
                        "Phy_Voxel_Straight", // its name
                        logicDNAEnvelope,     // its mother volume
                        false,                // no use
                        i,                    // copy number
                        fCheckOverlaps);      // checking overlaps
    AddStructureToMap(i, structure, fOffsets[i - 1]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  G4cout << "construction time: " << duration.count() << " ms" << G4endl;

  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField() {
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare World as a MultiFunctionalDetector scorer
  G4MultiFunctionalDetector *mfDetector =
      new G4MultiFunctionalDetector("mfDetector");

  // LET scorer:
  //  - scores restricted or unrestricted LET
  ScoreLET *LET = new ScoreLET("LET");
  mfDetector->RegisterPrimitive(LET);

  // Species scorer:
  //  - scores number of species over time
  //  - compute the radiochemical yields (G values)
  G4VPrimitiveScorer *primitivSpecies = new ScoreSpecies("Species");
  mfDetector->RegisterPrimitive(primitivSpecies);

  // SB Scorer: Requires access to Geometry
  //  - scores number of Direct SB
  //  - scores number of Indirect SB
  G4VPrimitiveScorer *primitiveSB =
      new ScoreStrandBreaks("StrandBreaks", this, fEnvelopeDiameter);
  mfDetector->RegisterPrimitive(primitiveSB);

  // Attach Detectors to DNA Volumes and register to SDManager
  G4LogicalVolumeStore *theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t t = 0; t < theLogicalVolumes->size(); t++) {
    G4String lVolumeName = (*theLogicalVolumes)[t]->GetName();
    // ! attach SB detector to all logical volumes except the world volume
    if (lVolumeName != "Log_World") {
      (*theLogicalVolumes)[t]->SetSensitiveDetector(mfDetector);
    }
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(mfDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ReadOffsetFile(G4String fileName) {
  if (fOffsets.size() > 0)
    fOffsets.clear();

  std::ifstream offsetFile;
  offsetFile.open(fileName);
  G4double x, y, z;

  if (!offsetFile) {
    G4cout << "DNA Offset positions file not found!!!" << G4endl;
    exit(1);
  } else {
    while (offsetFile >> x >> y >> z) {
      fOffsets.push_back(
          G4ThreeVector(x * angstrom, y * angstrom, z * angstrom));
    }
  }
  offsetFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::AddStructureToMap(G4int structureId,
                                             DNAStructure &structure,
                                             G4ThreeVector offset) {
  // Move DNA structure to the specified offset position and add to the map
  DNAStructure movedStructure =
      DNAHandler::GetMovedStructure(structureId, structure, offset);
  fStructureMap[structureId] = movedStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DNAStructure& DetectorConstruction::GetStructure(G4int structureId) {
  if (fStructureMap.find(structureId) == fStructureMap.end()) {
    G4cout << "DNA DNAStructure with ID " << structureId << " not found!"
           << G4endl;
    exit(1);
  }
  return fStructureMap[structureId];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
