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

#include "Run.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "G4VSensitiveDetector.hh"

#include "RunAction.hh"
#include "ScoreSpecies.hh"
#include "TimeStepAction.hh"

#include <iterator>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

Run::Run() : G4Run() {
  G4MultiFunctionalDetector *mfdet = dynamic_cast<G4MultiFunctionalDetector *>(
      G4SDManager::GetSDMpointer()->FindSensitiveDetector("mfDetector"));
  fCollectionIDSpecies =
      G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/Species");
  fCollectionIDLET =
      G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/LET");
  G4int CollectionSB =
      G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/StrandBreaks");

  fTotalLET = new G4THitsMap<G4double>("mfDetector", "LET");
  fScorerRun = mfdet->GetPrimitive(fCollectionIDSpecies);
  fLETScorerRun = mfdet->GetPrimitive(fCollectionIDLET);
  fStrandBreakRun = mfdet->GetPrimitive(CollectionSB);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Run::RecordEvent(const G4Event *event) {
  if (event->IsAborted())
    return;

  // Hits collections
  //
  G4HCofThisEvent *HCE = event->GetHCofThisEvent();
  if (!HCE)
    return;

  G4THitsMap<G4double> *evtMap =
      static_cast<G4THitsMap<G4double> *>(HCE->GetHC(fCollectionIDSpecies));

  G4THitsMap<G4double> *evtLET =
      static_cast<G4THitsMap<G4double> *>(HCE->GetHC(fCollectionIDLET));

  G4int nOfMap = evtLET->entries();

  G4int nOftotal = fTotalLET->entries();

  for (G4int i = 0; i < nOfMap; i++) {
    G4double *LET = (*evtLET)[i];
    if (!LET)
      continue;
    fTotalLET->add(nOftotal + i, *LET);
  }

  std::map<G4int, G4double *>::iterator itr;

  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep = *(itr->second);
    fSumEne += edep;
  }

  // Collect species records from TimeStepAction at end of each event
  if (G4Scheduler::Instance() != nullptr) {
    TimeStepAction *timeStepAction = dynamic_cast<TimeStepAction *>(
        G4Scheduler::Instance()->GetUserTimeStepAction());
    if (timeStepAction) {
      std::vector<SpeciesRecord> records = timeStepAction->GetSpeciesRecords();
      for (const auto &record : records) {
        fSpeciesRecords.push_back(record);
      }
      timeStepAction->ClearSpeciesRecords(); // Clear after collecting
    }
  }

  G4Run::RecordEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Run::Merge(const G4Run *aRun) {
  if (aRun == this) {
    return;
  }

  const Run *localRun = static_cast<const Run *>(aRun);
  fSumEne += localRun->fSumEne;

  G4int nOfMaster = fTotalLET->entries();
  G4int nOfLocal = localRun->fTotalLET->entries();
  for (G4int i = 0; i < nOfLocal; i++) {
    G4double *LET = (*localRun->fTotalLET)[i];
    if (!LET)
      continue;
    fTotalLET->add(nOfMaster + i, *LET);
  }

  ScoreSpecies *masterSpeciesScorer =
      dynamic_cast<ScoreSpecies *>(this->fScorerRun);

  ScoreSpecies *localSpeciesScorer =
      dynamic_cast<ScoreSpecies *>(localRun->fScorerRun);

  masterSpeciesScorer->AbsorbResultsFromWorkerScorer(localSpeciesScorer);

  ScoreStrandBreaks *masterSBScorer =
      dynamic_cast<ScoreStrandBreaks *>(this->fStrandBreakRun);

  ScoreStrandBreaks *localSBScorer =
      dynamic_cast<ScoreStrandBreaks *>(localRun->fStrandBreakRun);

  masterSBScorer->AbsorbResultsFromWorkerScorer(localSBScorer);

  for (const auto &record : localRun->fSpeciesRecords) {
    this->fSpeciesRecords.push_back(record);
  }

  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
