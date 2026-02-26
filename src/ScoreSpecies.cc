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

#include "ScoreSpecies.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Scheduler.hh"
#include "G4TScoreNtupleWriter.hh"
#include "G4UImessenger.hh"
#include "G4UnitsTable.hh"

#include <G4EventManager.hh>
#include <G4MolecularConfiguration.hh>
#include <G4MoleculeCounter.hh>
#include <G4SystemOfUnits.hh>
#include <globals.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::ScoreSpecies(G4String name, G4int depth)
  : G4VPrimitiveScorer(name, depth),
    G4UImessenger(),
    fOutputType("root"),  // other options: "csv", "hdf5", "xml"
    fHCID(-1)
{
  fpSpeciesdir = new G4UIdirectory("/scorer/species/");
  fpSpeciesdir->SetGuidance("ScoreSpecies commands");

  fpAddTimeToRecordcmd = new G4UIcmdWithADoubleAndUnit("/scorer/species/addTimeToRecord", this);

  fpTimeBincmd = new G4UIcmdWithAnInteger("/scorer/species/setNumberOfTimeBins", this);

  fpOutputTypeUI = new G4UIcmdWithAString("/scorer/species/setOutputFormat", this);
  fpOutputFileUI = new G4UIcmdWithAString("/scorer/species/setOutputFile", this);

  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::~ScoreSpecies()
{
  delete fpSpeciesdir;
  delete fpAddTimeToRecordcmd;
  delete fpTimeBincmd;
  delete fpOutputTypeUI;
  delete fpOutputFileUI;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpAddTimeToRecordcmd)
  {
    G4double cmdTime = fpAddTimeToRecordcmd->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if (command == fpTimeBincmd)
  {
    ClearTimeToRecord();
    G4int cmdBins = fpTimeBincmd->GetNewIntValue(newValue);
    G4double timeMin = 1 * ps;
    G4double timeMax = G4Scheduler::Instance()->GetEndTime() - 1 * ps;
    G4double timeLogMin = std::log10(timeMin);
    G4double timeLogMax = std::log10(timeMax);
    for (G4int i = 0; i < cmdBins; i++)
    {
      AddTimeToRecord(std::pow(10, timeLogMin + i * (timeLogMax - timeLogMin) / (cmdBins - 1)));
    }
  }

  if (command == fpOutputTypeUI)
  {
    fOutputType = newValue;
    G4StrUtil::to_lower(fOutputType);
  }

  if (command == fpOutputFileUI) fOutputFile = newValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreSpecies::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.) return FALSE;

  edep *= aStep->GetPreStepPoint()->GetWeight();
  G4int index = GetIndex(aStep);
  fEvtMap->add(index, edep);
  fEdep += edep;

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::Initialize(G4HCofThisEvent* HCE)
{
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());

  if (fHCID < 0)
  {
    fHCID = GetCollectionID(0);
  }

  HCE->AddHitsCollection(fHCID, (G4VHitsCollection*)fEvtMap);
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::EndOfEvent(G4HCofThisEvent*)
{
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted())
  {
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }

  G4MoleculeCounter::RecordedMolecules species =
    G4MoleculeCounter::Instance()->GetRecordedMolecules();

  if (species.get() == 0 || species->size() == 0)
  {
    G4cout << "No molecule recorded, energy deposited= " << G4BestUnit(fEdep, "Energy") << G4endl;
    ++fNEvent;
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }

  for (auto molecule : *species)
  {
    for (auto time_mol : fTimeToRecord)
    {
      double n_mol = G4MoleculeCounter::Instance()->GetNMoleculesAtTime(molecule, time_mol);

      if (n_mol < 0)
      {
        G4cerr << "N molecules not valid < 0 " << G4endl;
        G4Exception("", "N<0", FatalException, "");
      }
      // note: fSpeciesInfoPerTime[time_mol] will create the InnerSpeciesMap if it does not exist
      // then empty SpeciesInfo will be created for molecule (all members zero)
      SpeciesInfo& molInfo = fSpeciesInfoPerTime[time_mol][molecule];
      molInfo.fNumber += n_mol;
      molInfo.fNumber2 += n_mol * n_mol;

      double gValue = 0.0;
      if (fEdep > 0.)
      {
        gValue = (n_mol / (fEdep / eV)) * 100.;
      }

      molInfo.fG += gValue;
      molInfo.fG2 += gValue * gValue;
    }
  }

  ++fNEvent;

  fEdep = 0.;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreSpecies* right = dynamic_cast<ScoreSpecies*>(workerScorer);

  if (right == 0)
  {
    return;
  }
  if (right == this)
  {
    return;
  }

  SpeciesMap::iterator it_map1 = right->fSpeciesInfoPerTime.begin();
  SpeciesMap::iterator end_map1 = right->fSpeciesInfoPerTime.end();

  for (; it_map1 != end_map1; ++it_map1)
  {
    InnerSpeciesMap& map2 = it_map1->second;
    InnerSpeciesMap::iterator it_map2 = map2.begin();
    InnerSpeciesMap::iterator end_map2 = map2.end();

    for (; it_map2 != end_map2; ++it_map2)
    {
      SpeciesInfo& molInfo = fSpeciesInfoPerTime[it_map1->first][it_map2->first];
      molInfo.fNumber += it_map2->second.fNumber;
      molInfo.fNumber2 += it_map2->second.fNumber2;
      molInfo.fG += it_map2->second.fG;
      molInfo.fG2 += it_map2->second.fG2;
    }
  }
  right->fSpeciesInfoPerTime.clear();

  fNEvent += right->fNEvent;
  right->fNEvent = 0;
  right->fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::PrintAll()
{
  G4cout << " MultiFunctionalDet    " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of events " << fNEvent << G4endl;
  G4cout << " Number of energy deposition recorded " << fEvtMap->entries() << G4endl;

  for (auto itr : *fEvtMap->GetMap())
  {
    G4cout << "    copy no.: " << itr.first
           << "    energy deposit: " << *(itr.second) / GetUnitValue() << " [" << GetUnit() << "]"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::OutputAndClear()
{
  if (G4Threading::IsWorkerThread()) return;

  //---------------------------------------------------------------------------
  // Save results

  if (fOutputType != "ascii")
  {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType(fOutputType);

    if (analysisManager)
    {
      this->WriteWithAnalysisManager(analysisManager);
    }
  }
  else
    OutputToASCII();

  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::WriteWithAnalysisManager(G4VAnalysisManager* analysisManager)
{
  analysisManager->OpenFile(fOutputFile+".root");
  int fNtupleID = analysisManager->CreateNtuple("species", "species");

  // Match columns to OutputToASCII: Time, G, RMS, Nb, Molecule
  analysisManager->CreateNtupleDColumn(fNtupleID, "Time");  // Time [ps]
  analysisManager->CreateNtupleDColumn(fNtupleID, "G");  // Gx(/100 eV)
  analysisManager->CreateNtupleDColumn(fNtupleID, "RMS");  // RMS
  analysisManager->CreateNtupleDColumn(fNtupleID, "Nb");  // Gx(#Mols)
  analysisManager->CreateNtupleSColumn(fNtupleID, "Molecule");  // Molecule Name
  analysisManager->FinishNtuple(fNtupleID);

  for (auto it_map1 : fSpeciesInfoPerTime)
  {
    InnerSpeciesMap& map2 = it_map1.second;

    // Use ps for unit as in OutputToASCII
    double time = it_map1.first / ps;

    for (auto it_map2 : map2)
    {
      auto species = it_map2.first;
      const G4String& name = species->GetName();

      // Calculate statistical quantities
      double N = fNEvent;

      double sumNb = it_map2.second.fNumber;
      double sumNb2 =
        it_map2.second
          .fNumber2;
      double sumG = it_map2.second.fG;
      double sumG2 = it_map2.second.fG2;

      double avgNb = sumNb;
      double avgG = sumG;

      if (N > 0)
      {
        avgNb /= N;
        avgG /= N;
      }

      double sigmaG = sumG2;

      if (N == 1)
      {
        sigmaG = 0.0;
      }
      else if (N > 1)
      {
        sigmaG = std::sqrt(std::max(0.0, (N / (N - 1)) * ((sumG2 / N) - avgG * avgG)));
      }

      analysisManager->FillNtupleDColumn(fNtupleID, 0, time);
      analysisManager->FillNtupleDColumn(fNtupleID, 1, avgG);
      analysisManager->FillNtupleDColumn(fNtupleID, 2, sigmaG);
      analysisManager->FillNtupleDColumn(fNtupleID, 3, avgNb);
      analysisManager->FillNtupleSColumn(fNtupleID, 4, name);
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  analysisManager->Clear();
  fRunID++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::OutputToASCII()
{
  std::ofstream SpeciesOutput;
  SpeciesOutput.open(fOutputFile + ".txt");

  std::map<G4String, std::map<G4double, std::vector<G4double>>> mol;

  for (auto it_map1 : fSpeciesInfoPerTime)
  {
    InnerSpeciesMap& map2 = it_map1.second;
    G4double time = it_map1.first / ps;

    for (auto it_map2 : map2)
    {
      G4double G = it_map2.second.fG;
      G4double G2 = it_map2.second.fG2;
      G4double Nb = it_map2.second.fNumber;
      G4double Nb2 = it_map2.second.fNumber2;
      G4double N = fNEvent;
      if (N > 0)
      {
        G /= N;
        Nb /= N;
      }

      if (N == 1)
      {
        G2 = 0.0;
        Nb2 = 0.0;
      }

      else if (N > 0)
      {
        G2 = std::sqrt(N / (N - 1) * (G2 / N - G * G));
        Nb2 = std::sqrt(N / (N - 1) * (Nb2 / N - Nb * Nb));
      }

      mol[it_map2.first->GetName()][time].push_back(G);
      mol[it_map2.first->GetName()][time].push_back(G2);
      mol[it_map2.first->GetName()][time].push_back(Nb);
      mol[it_map2.first->GetName()][time].push_back(Nb2);
    }
  }

  SpeciesOutput << "# Species Scorer Output " << G4endl;
  SpeciesOutput << "# " << std::setw(10) << "Time [ps]" << std::setw(15) << "Gx(/100 eV)"
                << std::setw(15) << "RMS" << std::setw(15) << "Gx(#Mols)" << std::setw(20)
                << "Molecule" << G4endl;

  for (auto it1 : mol)
    for (auto it2 : it1.second)
      SpeciesOutput << std::setw(12) << it2.first << "     " << std::setw(12) << it2.second[0]
                    << "     " << std::setw(12) << it2.second[1] << "     " << std::setw(12)
                    << it2.second[2] << "     " << std::setw(17) << it1.first << G4endl;

  SpeciesOutput.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
