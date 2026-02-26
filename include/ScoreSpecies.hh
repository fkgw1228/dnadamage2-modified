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

#ifndef DNADAMAGE2_ScoreSpecies_h
#define DNADAMAGE2_ScoreSpecies_h 1

#include "G4THitsMap.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "G4VPrimitiveScorer.hh"

#include <set>

class G4VAnalysisManager;
class G4MolecularConfiguration;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ScoreSpecies : public G4VPrimitiveScorer, public G4UImessenger
{
  public:
    ScoreSpecies(G4String name, G4int depth = 0);
    ~ScoreSpecies() override;

    /** Add a time at which the number of species should be recorded.
        Default times are set up to 1 microsecond.*/
    inline void AddTimeToRecord(double time) { fTimeToRecord.insert(time); }

    /**  Remove all times to record, must be reset by user.*/
    inline void ClearTimeToRecord() { fTimeToRecord.clear(); }

    /** Get number of recorded events*/
    inline int GetNumberOfRecordedEvents() const { return fNEvent; }

    /** Write results to whatever chosen file format*/
    void WriteWithAnalysisManager(G4VAnalysisManager*);

    struct SpeciesInfo
    {
        SpeciesInfo() {}
        SpeciesInfo(const SpeciesInfo& right)  // Species A(B);
        {
          fNumber = right.fNumber;
          fNumber2 = right.fNumber2;
          fG = right.fG;
          fG2 = right.fG2;
        }
        SpeciesInfo& operator=(const SpeciesInfo& right)  // A = B
        {
          if (&right == this) return *this;
          fNumber2 = right.fNumber2;
          fNumber = right.fNumber;
          fG = right.fG;
          fG2 = right.fG2;
          return *this;
        }
        int fNumber = 0;
        int fNumber2 = 0;
        double fG = 0;
        double fG2 = 0;
    };

  private:
    typedef const G4MolecularConfiguration Species;
    typedef std::map<Species*, SpeciesInfo> InnerSpeciesMap;
    typedef std::map<double, InnerSpeciesMap> SpeciesMap;
    SpeciesMap fSpeciesInfoPerTime;

    std::set<G4double> fTimeToRecord;

    int fNEvent = 0;  // number of processed events
    double fEdep = 0;  // total energy deposition
    G4String fOutputType;  // output type

  protected:
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

  public:
    void Initialize(G4HCofThisEvent*) override;
    void EndOfEvent(G4HCofThisEvent*) override;
    void DrawAll() override;
    void PrintAll() override;
    /** Method used in multithreading mode in order to merge
        the results*/
    void AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer*);
    void OutputAndClear();
    void SetNewValue(G4UIcommand*, G4String) override;
    void OutputToASCII();

    SpeciesMap GetSpeciesInfo() { return fSpeciesInfoPerTime; }

  private:
    G4int fHCID;
    G4THitsMap<G4double>* fEvtMap = nullptr;

    G4int fRunID = 0;
    G4UIdirectory* fpSpeciesdir = nullptr;
    G4UIcmdWithAnInteger* fpTimeBincmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fpAddTimeToRecordcmd = nullptr;

    G4UIcmdWithAString* fpOutputTypeUI = nullptr;
    G4UIcmdWithAString* fpOutputFileUI = nullptr;
    G4String fOutputFile = "Species_Info";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
