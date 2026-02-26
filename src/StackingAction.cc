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

#include "StackingAction.hh"

#include "G4DNAChemistryManager.hh"
#include "G4EventManager.hh"
#include "G4ITTrackHolder.hh"
#include "G4ManyFastLists.hh"
#include "G4SDManager.hh"
#include "G4Scheduler.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "TimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

StackingAction::StackingAction() : G4UserStackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

StackingAction::~StackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void StackingAction::NewStage() {
  if (stackManager->GetNTotalTrack() == 0) {
    // Set species data collection flag based on event ID
    if (G4Scheduler::Instance() != nullptr) {
      TimeStepAction *timeStepAction = dynamic_cast<TimeStepAction *>(
          G4Scheduler::Instance()->GetUserTimeStepAction());
      if (timeStepAction) {
        G4int eventID = G4EventManager::GetEventManager()
                            ->GetConstCurrentEvent()
                            ->GetEventID();
        timeStepAction->SetCollectSpeciesData(
            eventID < RunAction::GetRecordEventNumToRecord());
      }
    }
    G4TrackManyList *trackList = G4ITTrackHolder::Instance()->GetMainList();
    G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
    G4ManyFastLists<G4Track>::iterator it_end = trackList->end();
    for (; it_begin != it_end; ++it_begin) {
      if (it_begin->GetGlobalTime() < 0.) {
        it_begin->SetGlobalTime(G4Scheduler::Instance()->GetGlobalTime());
      }
    }
    // Run Chemistry
    G4DNAChemistryManager::Instance()->Run();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
