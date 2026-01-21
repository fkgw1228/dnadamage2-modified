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
// This example is provided by the Geant4-DNA collaboration
// dnadamage2 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include <chrono>

#include "G4UImanager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "ScoreSpecies.hh"
#include "ScoreLET.hh"
#include "ScoreStrandBreaks.hh"

#include "DNAFileReader.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "DNAHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
    fpDetDir = new G4UIdirectory("/det/");
    fpDetDir->SetGuidance("Detector control.");

    fpWorldSizeUI = new G4UIcmdWithADoubleAndUnit("/det/WorldSize",this);
    fpWorldSizeUI->SetDefaultUnit("um");

    fpEnvelopeRadiusUI = new G4UIcmdWithADoubleAndUnit("/det/EnvelopeRadius",this);
    fpEnvelopeRadiusUI->SetDefaultUnit("um");

    fpOffSetFileUI = new G4UIcmdWithAString("/det/OffSetFile", this);
    fpDNANbUI      = new G4UIcmdWithAnInteger("/det/NbOfDNA", this);

    fpDNAFileUI   = new G4UIcmdWithAString("/det/DNAFile", this);
    fpUseDNAUI    = new G4UIcmdWithABool("/det/UseDNAVolumes", this);
    fpCheckOverlapsUI = new G4UIcmdWithABool("/det/CheckOverlaps", this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction()
{
    delete fpWorldSizeUI;
    delete fpEnvelopeRadiusUI;
    delete fpOffSetFileUI;
    delete fpDNAFileUI;
    delete fpUseDNAUI;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fpWorldSizeUI)
    {
        // Set the size of the world
        fWorldSize = fpWorldSizeUI->GetNewDoubleValue(newValue);
    }
    if(command == fpEnvelopeRadiusUI)
    {
        // Set the size of the envelope
        fEnvelopeRadius = fpEnvelopeRadiusUI->GetNewDoubleValue(newValue);
    }
    
    if (command == fpDNANbUI)
    {
        // Set the number of DNA molecules to be constructed
        fNbOfDNA = fpDNANbUI->GetNewIntValue(newValue);
    }

    if (command == fpOffSetFileUI) 
    {
        // Read the DNA offset file
        G4String fileName = newValue;
        ReadOffsetFile(fileName);
    }

    if (command == fpDNAFileUI)
    {
        // Set the DNA structure file name (json)
        fDNAFile = newValue;
    }

    if (command == fpUseDNAUI)
    {
        // Set whether to use DNA volumes or not
        fUseDNAVolumes = fpUseDNAUI->GetNewBoolValue(newValue);
    }

    if (command == fpCheckOverlapsUI)
    {
        // Set whether to check overlaps or not
        fCheckOverlaps = fpCheckOverlapsUI->GetNewBoolValue(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Water is defined from NIST material database
    G4NistManager * man = G4NistManager::Instance();
    G4Material* normalWater = man->FindOrBuildMaterial("G4_WATER");
    G4Material* worldWater  = man->BuildMaterialWithNewDensity("G4_WATER_WORLD",
                                                               "G4_WATER",
                                                               1.0*g/cm/cm/cm);

    // World
    G4Box* worldBox              = new G4Box("Sol_World",
                                             0.5*fWorldSize,
                                             0.5*fWorldSize,
                                             0.5*fWorldSize);
    G4LogicalVolume* logicWorld  = new G4LogicalVolume(worldBox, 
                                                       worldWater, 
                                                       "Log_World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0,               // no rotation
                                                     G4ThreeVector(), // at (0,0,0)
                                                     "Phy_World",        // name of physical volume
                                                     logicWorld,      // logical volume
                                                     0,               // no mother volume
                                                     false,           // no use
                                                     0,               // copy number 
                                                     false);          // checking overlaps     
    
    // World Vis Attributes
    G4VisAttributes* worldVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    worldVis->SetVisibility(false);
    logicWorld->SetVisAttributes(worldVis);

    // Envelope for DNA geometry
    G4Orb* dnaEnvelope = new G4Orb("Sol_DNA_Envelope", fEnvelopeRadius);
    G4LogicalVolume* logicDNAEnvelope = new G4LogicalVolume(dnaEnvelope, 
                                                            normalWater,
                                                            "Log_DNA_Envelope");
    new G4PVPlacement(0,                  // no rotation
                      G4ThreeVector(),    // at (0,0,0)
                      logicDNAEnvelope,   // its logical volume
                      "Phy_DNA_Envelope", // its name  
                      logicWorld,         // its mother  volume
                      false,              // no use
                      0,                  // copy number
                      false);             // checking overlaps

    // DNA Envelope Vis Attributes
    G4VisAttributes* dnaEnvelopeVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    dnaEnvelopeVis->SetForceSolid(true);
    logicDNAEnvelope->SetVisAttributes(dnaEnvelopeVis);

    // DNA geometry construction
    auto start = std::chrono::high_resolution_clock::now();    // timer start
    DNAFileReader reader = DNAFileReader(fDNAFile);
    DNAStructure dnaStructure = reader.ReadDNAFile(1);    // Construct DNA structure from json file
    DNAHandler::Print(dnaStructure);

    // Voxel to store DNA geometry
    G4LogicalVolume* logicStraightVoxel;
    DNAGeometryConstructor constructor = DNAGeometryConstructor();
    logicStraightVoxel = constructor.CreateDNAGeometry(dnaStructure, fCheckOverlaps);

    for (G4int i=1; i<=fNbOfDNA; i++) {
        if (fUseDNAVolumes)
            new G4PVPlacement(0,                      // no rotation
                              fOffsets[i-1],          // at offset position 
                              logicStraightVoxel,     // its logical volume
                              "Phy_Voxel_Straight",   // its name
                              logicDNAEnvelope,       // its mother volume
                              false,                  // no use
                              i,                      // copy number
                              false);                 // checking overlaps
        AddDNAToMap(i, dnaStructure, fOffsets[i-1]);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    G4cout << "construction time: " << duration.count() << " ms" << G4endl;

    //always return the physical World
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    // declare World as a MultiFunctionalDetector scorer
    G4MultiFunctionalDetector* mfDetector =
            new G4MultiFunctionalDetector("mfDetector");

    // LET scorer:
    //  - scores restricted or unrestricted LET
    ScoreLET* LET = new ScoreLET("LET");
    mfDetector->RegisterPrimitive(LET);

    // Species scorer:
    //  - scores number of species over time
    //  - compute the radiochemical yields (G values)
    G4VPrimitiveScorer* primitivSpecies = new ScoreSpecies("Species");
    mfDetector->RegisterPrimitive(primitivSpecies);

    // SB Scorer: Requires access to Geometry
    //  - scores number of Direct SB
    //  - scores number of Indirect SB
    G4VPrimitiveScorer* primitiveSB = new ScoreStrandBreaks("StrandBreaks",
                                                          this, &fEnvelopeRadius);
    mfDetector->RegisterPrimitive(primitiveSB);

    // Attach Detectors to DNA Volumes and register to SDManager
    G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
    for ( size_t t = 0; t < theLogicalVolumes->size(); t++ ) {
        G4String lVolumeName = (*theLogicalVolumes)[t]->GetName();
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
    }
    else {
        while (offsetFile >> x >> y >> z) {
            fOffsets.push_back(G4ThreeVector(x*angstrom,y*angstrom,z*angstrom));
        }
    }  
    offsetFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::AddDNAToMap(G4int dnaId, DNAStructure& dnaStructure, G4ThreeVector offset) {
    DNAStructure movedDnaStructure = DNAHandler::GetMovedDNA(dnaId, dnaStructure, offset);
    fDNAMap[dnaId] = movedDnaStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const DNAStructure& DetectorConstruction::GetDnaStructure(G4int dnaId) {
    if (fDNAMap.find(dnaId) == fDNAMap.end()) {
        G4cout << "DNA Structure with ID " << dnaId << " not found!" << G4endl;
        exit(1);
    }
    return fDNAMap[dnaId];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....