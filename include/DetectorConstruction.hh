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
// DNADAMAGE2 example is derived from the chem6 example
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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DNADAMAGE2_DetectorConstruction_h
#define DNADAMAGE2_DetectorConstruction_h 1

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "Structure.hh"

#include <map>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction, public G4UImessenger
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    void SetNewValue(G4UIcommand*, G4String) override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    Structure& GetStructure(G4int);
    std::map<G4int, Structure> GetStructureMap() { return fStructureMap; }

  private:
    void ReadOffsetFile(G4String);
    void AddStructureToMap(G4int, Structure&, G4ThreeVector);

    // UI commands
    G4UIdirectory* fpDetDir = nullptr;
    G4UIcmdWithAString* fpOffsetFileUI = nullptr;
    G4UIcmdWithAString* fpDNAFileUI = nullptr;
    G4UIcmdWithAnInteger* fpNumberOfDNAUI = nullptr;
    G4UIcmdWithADoubleAndUnit* fpWorldSizeUI = nullptr;
    G4UIcmdWithADoubleAndUnit* fpEnvelopeDiameterUI = nullptr;
    G4UIcmdWithABool* fpEnableDNAVolumesUI = nullptr;
    G4UIcmdWithABool* fpCheckOverlapsUI = nullptr;

    // Geometry parameters
    G4double fWorldSize = 1 * um;
    G4double fEnvelopeDiameter = 0.5 * um;

    std::vector<G4ThreeVector> fOffsets;      // Offsets for DNA placements
    G4String fDNAFile = "";                   // DNA structure file name
    G4int fNumberOfDNA = 0;                   // Number of DNA structures to be placed
    G4bool fEnableDNAVolumes = false;         // Flag to use DNA volumes for simulation
    G4bool fCheckOverlaps = false;            // Flag to check overlaps during geometry construction
    std::map<G4int, Structure> fStructureMap; // Placed DNA structures keyed by their IDs
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
