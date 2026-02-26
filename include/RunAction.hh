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

#ifndef DNADAMAGE2_RunAction_h
#define DNADAMAGE2_RunAction_h 1

#include "G4THitsMap.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "G4UserRunAction.hh"

#include "TimeStepAction.hh"

class G4Run;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction, public G4UImessenger {
public:
  RunAction();
  // TIPs: please avoid constructors with arguments
  // all data can be retrieved from G4RunManager
  // or others: G4SDManager::FindSensitiveDetector
  ~RunAction() override;

  void SetNewValue(G4UIcommand *, G4String) override;
  static G4int GetRecordEventNumToRecord() { return fgRecordEventNumToRecord; }

  void OutputSpeciesPositions(const std::vector<SpeciesRecord> &speciesRecords);

  G4Run *GenerateRun() override;
  void BeginOfRunAction(const G4Run *) override;
  void EndOfRunAction(const G4Run *) override;

private:
  G4UIdirectory *fpSpeciesPosDir = nullptr;
  G4UIcmdWithAString *fpOutputFileNameUI = nullptr;
  G4UIcmdWithAString *fpOutputTypeUI = nullptr;
  G4UIcmdWithAnInteger *fpRecordNumUI = nullptr;

  G4String fOutputFileName = "species_positions";
  G4String fOutputType = "root";
  static G4int fgRecordEventNumToRecord;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
