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

#ifndef DNADAMAGE2_ScoreLET_h
#define DNADAMAGE2_ScoreLET_h 1

#include "G4THitsMap.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "G4VPrimitiveScorer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ScoreLET : public G4VPrimitiveScorer, public G4UImessenger
{
  public:  // with description
    ScoreLET(G4String name);
    ~ScoreLET() override;
    void Initialize(G4HCofThisEvent*) override;
    void EndOfEvent(G4HCofThisEvent*) override;
    void OutputAndClear();

    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    G4int GetIndex(G4Step*) override;
    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    G4UIdirectory* fpLETDir = nullptr;
    G4UIcmdWithADoubleAndUnit* fpCutoff = nullptr;

    G4double fCutoff = DBL_MAX;
    G4int fNEvent = 0;
    G4double fLET = 0;
    G4double fEdep = 0;
    G4double fStepL = 0;
    G4int fTrackID = 1;

    G4THitsMap<G4double>* fEvtMap = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
