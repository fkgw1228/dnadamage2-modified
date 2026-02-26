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

#ifndef DNADAMAGE2_TimeStepAction_h
#define DNADAMAGE2_TimeStepAction_h 1

#include "G4SystemOfUnits.hh"
#include "G4UserTimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct SpeciesRecord
{
    G4int eventID;
    G4int trackID;
    G4String speciesName;
    G4double time;
    G4ThreeVector position;

    bool operator==(const SpeciesRecord& other) const
    {
      return (eventID == other.eventID && trackID == other.trackID
              && speciesName == other.speciesName && time == other.time);
    }
};

typedef std::map<G4int, SpeciesRecord> SpeciesRecordPerTrackMap;
typedef std::map<G4int, SpeciesRecordPerTrackMap> SpeciesRecordPerEventMap;

class TimeStepAction : public G4UserTimeStepAction
{
  public:
    TimeStepAction();
    ~TimeStepAction() override = default;
    TimeStepAction(const TimeStepAction& other);
    TimeStepAction& operator=(const TimeStepAction& other);

    std::vector<SpeciesRecord> GetSpeciesRecords() const { return fSpeciesRecords; }
    void ClearSpeciesRecords() { fSpeciesRecords.clear(); }
    void SetCollectSpeciesData(G4bool collect) { fCollectSpeciesData = collect; }
    G4bool GetCollectSpeciesData() const { return fCollectSpeciesData; }

    void StartProcessing() override {}
    void UserPreTimeStepAction() override;
    void UserPostTimeStepAction() override;

    void UserReactionAction(const G4Track& /*trackA*/, const G4Track& /*trackB*/,
                            const std::vector<G4Track*>* /*products*/) override;

    void EndProcessing() override {}
    void Clear() {}

  private:
    std::vector<SpeciesRecord> fSpeciesRecords;
    G4bool fCollectSpeciesData = false;
    const G4double kXRange = 100 * angstrom;
    const G4double kYRange = 100 * angstrom;
    const G4double kZRange = 150 * angstrom;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
