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

#ifndef DNADAMAGE2_PhysicsList_h
#define DNADAMAGE2_PhysicsList_h 1

#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserChemistryList.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include <G4UImessenger.hh>

class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList : public G4VModularPhysicsList, public G4UImessenger
{
  public:
    explicit PhysicsList();
    ~PhysicsList() override;

    void ConstructParticle() override;
    void ConstructProcess() override;

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

    void RegisterPhysicsConstructor(const G4String& name);
    void RegisterChemistryConstructor(const G4String& name);

  private:
    G4VPhysicsConstructor* fEmDNAPhysicsList = nullptr;
    G4VPhysicsConstructor* fEmDNAChemistryList = nullptr;
    G4String fPhysDNAName = "";
    G4String fChemDNAName = "";

    // UI commands
    G4UIcmdWithAString* fpPhysicsUI = nullptr;
    G4UIcmdWithAString* fpChemistryUI = nullptr;

    G4UIcmdWithADouble* fpDMSOUI = nullptr;
    G4UIcmdWithADouble* fpOxygenUI = nullptr;

    G4double fDMSO = 1E-2;
    G4double fOxygen = 0.27E-3;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
