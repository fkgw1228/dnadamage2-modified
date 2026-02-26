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

#include "ITSteppingAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITSteppingAction::ITSteppingAction() : G4UserSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITSteppingAction::~ITSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* track = aStep->GetTrack();
  const G4String& moleculeName = track->GetDefinition()->GetParticleName();

  G4bool isDNAMolecule = G4StrUtil::contains(moleculeName, "DNA_Deoxyribose")
                         || G4StrUtil::contains(moleculeName, "DNA_Histone");
  if (isDNAMolecule)
  {
    return;
  }

  const G4String& volumeName = preStepPoint->GetPhysicalVolume()->GetName();
  G4bool isInHistone = false;
  for (const G4String& aminoAcidName : fAminoAcidNames)
  {
    if (G4StrUtil::contains(volumeName, aminoAcidName))
    {
      isInHistone = true;
      break;
    }
  }
  // remove non-DNA molecules inside DNA and histone volumes
  if (isInHistone && !isDNAMolecule)
  {
    track->SetTrackStatus(fStopAndKill);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
