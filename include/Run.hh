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

#ifndef DNADAMAGE2_Run_h
#define DNADAMAGE2_Run_h 1

#include "G4Run.hh"
#include "G4THitsMap.hh"

#include "ScoreSpecies.hh"
#include "ScoreStrandBreaks.hh"
#include "TimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VPrimitiveScorer;
class Run : public G4Run {
public:
  Run();
  ~Run() override = default;

  void RecordEvent(const G4Event *) override;
  void Merge(const G4Run *) override;

  G4double GetSumDose() const { return fSumEne; }
  G4VPrimitiveScorer *GetPrimitiveScorer() const { return fScorerRun; }
  G4VPrimitiveScorer *GetSBScorer() const { return fStrandBreakRun; }
  G4THitsMap<G4double> *GetLET() { return fTotalLET; }

  const std::vector<SpeciesRecord> &GetSpeciesRecords() const {
    return fSpeciesRecords;
  }
  void AddSpeciesRecord(SpeciesRecord record) {
    fSpeciesRecords.push_back(record);
  }

private:
  G4double fSumEne = 0;
  G4VPrimitiveScorer *fScorerRun = nullptr;
  G4VPrimitiveScorer *fLETScorerRun = nullptr;
  G4VPrimitiveScorer *fStrandBreakRun = nullptr;
  G4THitsMap<G4double> *fTotalLET = nullptr;
  G4int fCollectionIDSpecies = -1;
  G4int fCollectionIDLET = -1;

  std::vector<SpeciesRecord> fSpeciesRecords;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
