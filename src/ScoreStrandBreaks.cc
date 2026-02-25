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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file ScoreStrandBreaks.cc
/// \brief Implementation of the ScoreStrandBreaks class

#include "ScoreStrandBreaks.hh"

#include "G4AnalysisManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Scheduler.hh"
#include "G4UImessenger.hh"
#include "G4UnitsTable.hh"

#include "Model.hh"
#include "MoleculeInserter.hh"
#include "DNAStructure.hh"
#include "TimeStepAction.hh"
#include <G4EventManager.hh>
#include <G4SystemOfUnits.hh>
#include <globals.hh>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreStrandBreaks::ScoreStrandBreaks(G4String name, DetectorConstruction* detect,
                                     G4double envelopeDiameter)
  : G4VPrimitiveScorer(name), G4UImessenger(), fpDetector(detect)
{
  fpOutputFileUI = new G4UIcmdWithAString("/scorer/StrandBreak/setOutputFile", this);
  fpOutputFileUI->SetGuidance("Set output file name");

  fpOutputTypeUI = new G4UIcmdWithAString("/scorer/StrandBreak/setOutputFormat", this);
  fpOutputTypeUI->SetGuidance("Set output file format: ASCII by default");

  fpBreakEnergyUI = new G4UIcmdWithADoubleAndUnit("/scorer/StrandBreak/setBreakEnergy", this);
  fpBreakEnergyUI->SetDefaultUnit("eV");
  fpBreakEnergyUI->SetGuidance("Direct SSB energy threshold");

  fEnvelopeDiameter = envelopeDiameter;
  fpGun = new MoleculeInserter(true);
  // G4DNAChemistryManager::Instance()->SetGun(fpGun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreStrandBreaks::~ScoreStrandBreaks()
{
  delete fpOutputFileUI;
  delete fpBreakEnergyUI;
  delete fpOutputTypeUI;
  delete fpGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpOutputFileUI)
  {
    fOutputName = newValue;
  }

  if (command == fpBreakEnergyUI)
  {
    fBreakEnergy = fpBreakEnergyUI->GetNewDoubleValue(newValue);
  }

  if (command == fpOutputTypeUI)
  {
    fOutputType = newValue;
    G4StrUtil::to_lower(fOutputType);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreStrandBreaks::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep <= 0)
  {
    return FALSE;
  }
  fEnergyDeposit += edep;

  // identify DNA component
  G4StepPoint* aStepPoint = aStep->GetTrack()->GetStep()->GetPreStepPoint();
  G4TouchableHandle aTouchable = aStepPoint->GetTouchableHandle();
  G4String volumeName = aTouchable->GetVolume()->GetName();

  G4bool isSugar = G4StrUtil::contains(volumeName, "deoxyribose");
  // G4bool isBase = G4StrUtil::contains(volumeName, "base");
  G4bool isPhosphate = G4StrUtil::contains(volumeName, "phosphate");

  // direct damage record
  if (isSugar || isPhosphate)
  {
    if (edep < fBreakEnergy)
    {
      return TRUE;
    }

    // volumeName format: "phy_ellipsoid_<chainId>_<residueId>_<compoundName>"
    // After SplitByUnderscore: [0]="phy", [1]="ellipsoid", [2]=chainId, [3]=residueId,
    // [4]=compoundName
    G4int eventId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    std::vector<G4String> tokens = SplitByUnderscore(volumeName);
    G4int modelId = aTouchable->GetReplicaNumber(1);
    G4String chainId = tokens[2];
    G4int residueId = std::stoi(tokens[3]);
    G4String compoundName = tokens[4];

    G4ThreeVector pos = aStepPoint->GetPosition();

    DNARecordKey location{eventId, modelId, chainId, residueId, compoundName};
    if (fDirectDamageMap.find(location) == fDirectDamageMap.end())
      fDirectDamageMap[location].push_back(
        GetApproximateAtomName(modelId, chainId, residueId, compoundName, pos));

    return TRUE;
  }

  // indirect damage preparation
  if (!fDNAInserted)
  {
    // Insert DNA molecules
    fpGun->Clean();
    fInsertedMolecules.clear();

    G4int eventId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    std::vector<G4String> aminoAcidNames = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                            "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                            "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

    // place H atoms in deoxyribose molecule
    for (auto& [structureId, structure] : fpDetector->GetStructureMap())
    {
      Model model = structure.GetModel();
      for (auto& [chainId, chain] : model.GetChains())
      {
        for (auto& [residueId, residue] : chain.GetResidues())
        {
          for (auto& [compoundName, compound] : residue.GetCompounds())
          {
            std::vector<Atom> atoms = compound.GetAtoms();
            DNARecordKey location{eventId, structureId, chainId, residueId, compoundName};

            if (compoundName == "deoxyribose")
            {
              for (Atom atom : atoms)
              {
                G4ThreeVector position = atom.GetPosition();
                G4String atom_name = atom.GetName();
                if (atom_name == "H1\'")
                  InsertMolecule("H1", location, position);
                else if (atom_name == "H2\'")
                  InsertMolecule("H2a", location, position);
                else if (atom_name == "H2\'\'")
                  InsertMolecule("H2b", location, position);
                else if (atom_name == "H3\'")
                  InsertMolecule("H3", location, position);
                else if (atom_name == "H4\'")
                  InsertMolecule("H4", location, position);
                else if (atom_name == "H5\'")
                  InsertMolecule("H5a", location, position);
                else if (atom_name == "H5\'\'")
                  InsertMolecule("H5b", location, position);
              }
            }
            else
            {
              if (eventId == 0)
              {
                if (std::find(aminoAcidNames.begin(), aminoAcidNames.end(), compoundName)
                    != aminoAcidNames.end())
                {
                  Ellipsoid ellipsoid = structure.GetEllipsoid(chainId, residueId, compoundName);
                  if (!ellipsoid.IsEmpty())
                  {
                    G4ThreeVector position = ellipsoid.GetCenter();
                    InsertMolecule("Histone", location, position);
                  }
                }
              }
            }
          }
        }
      }
    }
    fpGun->DefineTracks();
    fDNAInserted = true;
  }
  return FALSE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<G4String> ScoreStrandBreaks::SplitByUnderscore(G4String& s)
{
  std::vector<G4String> result;
  std::stringstream ss(s);
  G4String item;

  while (std::getline(ss, item, '_'))
  {
    result.push_back(item);
  }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4String ScoreStrandBreaks::GetApproximateAtomName(G4int modelId,
                                                   G4String chainId,
                                                   G4int residueId,
                                                   G4String compound_name,
                                                   G4ThreeVector position)
{
  DNAStructure structure = fpDetector->GetStructure(modelId);
  Model model = structure.GetModel();
  std::vector<Atom> atoms =
    model.GetChain(chainId).GetResidue(residueId).GetCompound(compound_name).GetAtoms();
  G4String closestAtomName = "";
  G4double minDistance = std::numeric_limits<G4double>::max();

  for (auto atom : atoms)
  {
    G4double radius = atom.GetRadius();
    G4double distanceToAtomSurface = (atom.GetPosition() - position).mag() - radius;
    if (distanceToAtomSurface < minDistance)
    {
      minDistance = distanceToAtomSurface;
      closestAtomName = atom.GetName();
    }
  }

  return closestAtomName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreStrandBreaks::InsertMolecule(G4String moleculeName, DNARecordKey location,
                                       G4ThreeVector position)
{
  fpGun->AddMolecule(moleculeName, position, 1.0 * ps);
  fInsertedMolecules.push_back(std::make_pair(location, moleculeName));
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreStrandBreaks::Initialize(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::EndOfEvent(G4HCofThisEvent*)
{
  // Get eventId
  G4int eventId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // Absorbed Dose
  G4double volume = 4.0 / 3 * 3.14159 * (fEnvelopeDiameter * 0.5) * (fEnvelopeDiameter * 0.5)
                    * (fEnvelopeDiameter * 0.5);
  G4double density = 1.0 * g / cm3;
  G4double mass = density * volume;
  G4double Dose = fEnergyDeposit / mass;
  fDoseArray[eventId] = Dose / gray;

  G4int direct_num = 0;
  G4int indirect_num = 0;

  // Indirect Strand Breaks
  std::vector<G4Track*> InsertedTracks = fpGun->GetInsertedTracks();

  for (size_t i = 0; i < InsertedTracks.size(); i++)
  {
    if (InsertedTracks[i]->GetTrackStatus() == fStopAndKill)
    {
      std::pair<DNARecordKey, G4String> insertedPositions = fInsertedMolecules[i];
      DNARecordKey location = insertedPositions.first;
      G4String moleculeName = insertedPositions.second;

      fIndirectDamageMap[location].push_back(moleculeName);
      indirect_num++;
    }
  }

  fEnergyDeposit = 0;
  fDNAInserted = false;
  fNbOfEvents++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreStrandBreaks* right =
    dynamic_cast<ScoreStrandBreaks*>(dynamic_cast<G4VPrimitiveScorer*>(workerScorer));

  if (right == 0) return;
  if (right == this) return;

  DNADamageMap WorkerDirectDamageMap = right->fDirectDamageMap;
  DNADamageMap WorkerIndirectDamageMap = right->fIndirectDamageMap;
  std::map<G4int, G4double> WorkerDose = right->fDoseArray;

  for (auto& damageLocation : WorkerDirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto atomName : damageLocation.second)
      fDirectDamageMap[location].push_back(atomName);
  }

  for (auto& damageLocation : WorkerIndirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto atomName : damageLocation.second)
      fIndirectDamageMap[location].push_back(atomName);
  }

  for (auto EventAndDose : WorkerDose)
  {
    G4int eventId = EventAndDose.first;
    G4double dose = EventAndDose.second;
    fDoseArray[eventId] = dose;
  }

  fNbOfEvents += right->fNbOfEvents;
  right->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::Clear()
{
  if (fDirectDamageMap.size() != 0) fDirectDamageMap.clear();

  if (fIndirectDamageMap.size() != 0) fIndirectDamageMap.clear();

  if (fDoseArray.size() != 0) fDoseArray.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::PrintAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::OutputAndClear(G4double LET, G4double LET_STD)
{
  if (G4Threading::IsWorkerThread()) return;

  if (fOutputType != "ascii")
  {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if (fOutputType == "csv")
      analysisManager->SetDefaultFileType("csv");
    else if (fOutputType == "root")
      analysisManager->SetDefaultFileType("root");
    else if (fOutputType == "xml")
      analysisManager->SetDefaultFileType("xml");

    WriteWithAnalysisManager(analysisManager, LET, LET_STD);
  }
  else
    ASCII(LET, LET_STD);

  Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::WriteWithAnalysisManager(G4VAnalysisManager* analysisManager, G4double LET,
                                                 G4double LET_STD)
{
  analysisManager->OpenFile(fOutputName+".root");
  int fNtupleID = analysisManager->CreateNtuple("Damage", "Direct And Indirect Damage");
  analysisManager->CreateNtupleIColumn(fNtupleID, "EventId");
  analysisManager->CreateNtupleIColumn(fNtupleID, "ModelId");
  analysisManager->CreateNtupleSColumn(fNtupleID, "ChainId");
  analysisManager->CreateNtupleIColumn(fNtupleID, "ResidueId");
  analysisManager->CreateNtupleSColumn(fNtupleID, "CompoundName");
  analysisManager->CreateNtupleSColumn(fNtupleID, "AtomOrMoleculeName");
  analysisManager->CreateNtupleIColumn(fNtupleID, "DamageType");
  analysisManager->FinishNtuple(fNtupleID);

  int fLetNtupleID = analysisManager->CreateNtuple("LET", "LET Information");
  analysisManager->CreateNtupleDColumn(fLetNtupleID, "meanLET");
  analysisManager->CreateNtupleDColumn(fLetNtupleID, "stdLET");
  analysisManager->FinishNtuple(fLetNtupleID);

  for (auto& damageLocation : fDirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto atomName : damageLocation.second)
    {
      analysisManager->FillNtupleIColumn(fNtupleID, 0, location.eventID);
      analysisManager->FillNtupleIColumn(fNtupleID, 1, location.modelID);
      analysisManager->FillNtupleSColumn(fNtupleID, 2, location.chainID);
      analysisManager->FillNtupleIColumn(fNtupleID, 3, location.residueID);
      analysisManager->FillNtupleSColumn(fNtupleID, 4, location.compoundName);
      analysisManager->FillNtupleSColumn(fNtupleID, 5, atomName);
      analysisManager->FillNtupleIColumn(fNtupleID, 6, 1);  // Direct Damage
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  for (auto& damageLocation : fIndirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto moleculeName : damageLocation.second)
    {
      analysisManager->FillNtupleIColumn(fNtupleID, 0, location.eventID);
      analysisManager->FillNtupleIColumn(fNtupleID, 1, location.modelID);
      analysisManager->FillNtupleSColumn(fNtupleID, 2, location.chainID);
      analysisManager->FillNtupleIColumn(fNtupleID, 3, location.residueID);
      analysisManager->FillNtupleSColumn(fNtupleID, 4, location.compoundName);
      analysisManager->FillNtupleSColumn(fNtupleID, 5, moleculeName);
      analysisManager->FillNtupleIColumn(fNtupleID, 6, 2);  // Indirect Damage
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->FillNtupleDColumn(fLetNtupleID, 0, LET);
  analysisManager->FillNtupleDColumn(fLetNtupleID, 1, LET_STD);
  analysisManager->AddNtupleRow(fLetNtupleID);

  analysisManager->Write();
  analysisManager->CloseFile();
  analysisManager->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreStrandBreaks::ASCII(G4double LET, G4double LET_STD)
{
  std::ofstream dna(fOutputName + ".txt");

  dna << "# LET = " << LET / (keV / um) << " +- " << LET_STD / (keV / um) << " keV / um" << G4endl;
  dna << "#" << std::setw(14) << "Event-ID" << std::setw(14) << "DNA-ID" << std::setw(14)
      << "Chain-ID" << std::setw(14) << "Residue-ID" << std::setw(20) << "Compound" << std::setw(20)
      << "Atom/Molecule" << std::setw(14) << "Break-ID" << std::setw(14) << "Dose (Gy) " << G4endl;

  for (auto damageLocation : fDirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto atomName : damageLocation.second)
    {
      dna << std::setw(14) << location.eventID << std::setw(14) << location.modelID << std::setw(14)
          << location.chainID << std::setw(14) << location.residueID << std::setw(20)
          << location.compoundName << std::setw(20) << atomName << std::setw(14) << 1
          << std::setw(14) << fDoseArray[location.eventID] << G4endl;
    }
  }

  for (auto damageLocation : fIndirectDamageMap)
  {
    DNARecordKey location = damageLocation.first;
    for (auto moleculeName : damageLocation.second)
    {
      dna << std::setw(14) << location.eventID << std::setw(14) << location.modelID << std::setw(14)
          << location.chainID << std::setw(14) << location.residueID << std::setw(20)
          << location.compoundName << std::setw(20) << moleculeName << std::setw(14) << 2
          << std::setw(14) << fDoseArray[location.eventID] << G4endl;
    }
  }

  dna.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......