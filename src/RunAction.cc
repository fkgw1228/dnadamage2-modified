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
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4MoleculeCounter.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "Run.hh"
#include "TimeStepAction.hh"

#include <algorithm>
#include <cctype>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4int RunAction::fgRecordEventNumToRecord = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::RunAction() : G4UserRunAction() {
  // Create UI directory
  fpSpeciesPosDir = new G4UIdirectory("/speciesPos/");
  fpSpeciesPosDir->SetGuidance("Species position output control");

  // Create UI commands
  fpOutputFileNameUI =
      new G4UIcmdWithAString("/speciesPos/setOutputFile", this);
  fpOutputFileNameUI->SetGuidance(
      "Set the output file name for species position data.");

  fpOutputTypeUI =
      new G4UIcmdWithAString("/speciesPos/setOutputFormat", this);
  fpOutputTypeUI->SetGuidance(
      "Set species position output format (ascii or root).");
  fpOutputTypeUI->SetCandidates("ascii root");

  fpRecordNumUI =
      new G4UIcmdWithAnInteger("/speciesPos/setEventNumToRecord", this);
  fpRecordNumUI->SetGuidance(
      "Set the number of events to record species position data.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::~RunAction() {
  delete fpRecordNumUI;
  delete fpOutputTypeUI;
  delete fpOutputFileNameUI;
  delete fpSpeciesPosDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fpOutputFileNameUI)
  {
    fOutputFileName = newValue;
  }
  else if (command == fpOutputTypeUI)
  {
    fOutputType = newValue;
    G4StrUtil::to_lower(fOutputType);
  }
  else if (command == fpRecordNumUI)
  {
    fgRecordEventNumToRecord = fpRecordNumUI->GetNewIntValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Run *RunAction::GenerateRun() {
  Run *run = new Run();
  return run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::BeginOfRunAction(const G4Run *) {
  // informs the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  Run *dnadamage3Run = static_cast<Run *>(const_cast<G4Run *>(run));

  // Species records are now collected in Run::RecordEvent after each event
  // This ensures data is available before Run::Merge is called

  // results
  //
  G4double sumDose = dnadamage3Run->GetSumDose();

  // print
  //
  if (IsMaster()) {
    G4cout << G4endl
           << "--------------------End of Global Run-----------------------"
           << G4endl << "  The run has " << nofEvents << " events " << G4endl;

    ScoreSpecies *masterScorer =
        dynamic_cast<ScoreSpecies *>(dnadamage3Run->GetPrimitiveScorer());

    ScoreStrandBreaks *masterSBScorer =
        dynamic_cast<ScoreStrandBreaks *>(dnadamage3Run->GetSBScorer());

    G4cout << "Number of events recorded by the species scorer = "
           << masterScorer->GetNumberOfRecordedEvents() << G4endl;

    // LET
    Run *aRun = (Run *)run;
    G4THitsMap<G4double> *totLET = aRun->GetLET();
    G4int nOfEvent = totLET->entries();
    G4double LET_mean = 0;
    G4double LET_square = 0;
    for (G4int i = 0; i < nOfEvent; i++) {
      G4double *LET = (*totLET)[i];
      if (!LET)
        continue;
      LET_mean += *LET;
      LET_square += (*LET) * (*LET);
    }
    LET_mean /= nOfEvent;
    LET_square = std::sqrt(LET_square / nOfEvent - std::pow(LET_mean, 2));

    masterScorer->OutputAndClear();
    masterSBScorer->OutputAndClear(LET_mean, LET_square);

    const std::vector<SpeciesRecord> &speciesRecords =
        dnadamage3Run->GetSpeciesRecords();
    G4cout << "Total species records to output: " << speciesRecords.size()
           << G4endl;
    OutputSpeciesPositions(speciesRecords);
  } else {
    G4cout << G4endl
           << "--------------------End of Local Run------------------------"
           << G4endl << "  The run has " << nofEvents << " events" << G4endl;
  }

  G4cout << " Total energy deposited in the world volume : " << sumDose / eV
         << " eV" << G4endl
         << " ------------------------------------------------------------"
         << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::OutputSpeciesPositions(
    const std::vector<SpeciesRecord> &speciesRecords) {
  if (speciesRecords.empty())
    return;

  if (fOutputType == "ascii") {
    std::ofstream outputFile(fOutputFileName + ".txt");
    if (!outputFile.is_open()) {
      G4cerr << "Failed to open species position output file: "
             << fOutputFileName << ".txt" << G4endl;
      return;
    }

    outputFile << "# EventID TrackID SpeciesName Time[ps] PosX[Ang] PosY[Ang]"
               << " PosZ[Ang]" << std::endl;

    for (const auto &record : speciesRecords) {
      outputFile << record.eventID << " " << record.trackID << " "
                 << record.speciesName << " " << record.time << " "
                 << record.position.x() << " " << record.position.y() << " "
                 << record.position.z() << std::endl;
    }
    outputFile.close();
    return;
  }

  // Create analysis manager
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->OpenFile(fOutputFileName+".root");
  G4int ntupleId = analysisManager->CreateNtuple(
      "SpeciesTimeDistribution", "Temporal Distribution of Chemical Species");
  analysisManager->CreateNtupleIColumn(ntupleId, "EventID");
  analysisManager->CreateNtupleIColumn(ntupleId, "TrackID");
  analysisManager->CreateNtupleSColumn(ntupleId, "SpeciesName");
  analysisManager->CreateNtupleDColumn(ntupleId, "Time");
  analysisManager->CreateNtupleDColumn(ntupleId, "PosX");
  analysisManager->CreateNtupleDColumn(ntupleId, "PosY");
  analysisManager->CreateNtupleDColumn(ntupleId, "PosZ");
  analysisManager->FinishNtuple(ntupleId);

  for (const auto &record : speciesRecords) {
    analysisManager->FillNtupleIColumn(ntupleId, 0, record.eventID);
    analysisManager->FillNtupleIColumn(ntupleId, 1, record.trackID);
    analysisManager->FillNtupleSColumn(ntupleId, 2, record.speciesName);
    analysisManager->FillNtupleDColumn(ntupleId, 3, record.time);
    analysisManager->FillNtupleDColumn(ntupleId, 4, record.position.x());
    analysisManager->FillNtupleDColumn(ntupleId, 5, record.position.y());
    analysisManager->FillNtupleDColumn(ntupleId, 6, record.position.z());
    analysisManager->AddNtupleRow(ntupleId);
  }
  analysisManager->Write();
  analysisManager->CloseFile();
  analysisManager->Clear();

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
