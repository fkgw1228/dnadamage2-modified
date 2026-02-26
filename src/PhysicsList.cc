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

#include "PhysicsList.hh"

#include "G4EmDNAChemistry.hh"
#include "G4EmDNAChemistry_option1.hh"
#include "G4EmDNAChemistry_option2.hh"
#include "G4EmDNAChemistry_option3.hh"
#include "G4EmDNAChemistry_userOption.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4SystemOfUnits.hh"

#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4double currentDefaultCut = 1. * nanometer;
  // fixe lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV, 1 * GeV);
  SetDefaultCutValue(currentDefaultCut);
  SetVerboseLevel(1);

  RegisterPhysicsConstructor("G4EmDNAPhysics_option2");
  RegisterChemistryConstructor("G4EmDNAChemistry_userOption");

  fpPhysicsUI = new G4UIcmdWithAString("/physics/setPhysics", this);
  fpPhysicsUI->AvailableForStates(G4State_PreInit);

  fpChemistryUI = new G4UIcmdWithAString("/physics/setChemistry", this);
  fpChemistryUI->AvailableForStates(G4State_PreInit);

  fpDMSOUI = new G4UIcmdWithADouble("/chem/scavenger/setDMSOMole", this);
  fpDMSOUI->AvailableForStates(G4State_PreInit);

  fpOxygenUI = new G4UIcmdWithADouble("/chem/scavenger/setOxygenMole", this);
  fpOxygenUI->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fpPhysicsUI;
  delete fpChemistryUI;
  delete fpDMSOUI;
  delete fpOxygenUI;
  delete fEmDNAPhysicsList;
  delete fEmDNAChemistryList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if (fEmDNAPhysicsList != nullptr)
  {
    fEmDNAPhysicsList->ConstructParticle();
  }
  if (fEmDNAChemistryList != nullptr)
  {
    fEmDNAChemistryList->ConstructParticle();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  if (fEmDNAPhysicsList != nullptr)
  {
    fEmDNAPhysicsList->ConstructProcess();
  }
  if (fEmDNAChemistryList != nullptr)
  {
    fEmDNAChemistryList->ConstructProcess();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterPhysicsConstructor(const G4String& name)
{
  if (name == fPhysDNAName)
  {
    return;
  }

  if (name == "G4EmDNAPhysics")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option1")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option1(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option2")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option2(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option3")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option3(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option4")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option4(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option5")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option5(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option6")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option6(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option7")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option7(verboseLevel);
    fPhysDNAName = name;
  }
  else if (name == "G4EmDNAPhysics_option8")
  {
    delete fEmDNAPhysicsList;
    fEmDNAPhysicsList = new G4EmDNAPhysics_option8(verboseLevel);
    fPhysDNAName = name;
  }
  else
  {
    G4cout << "PhysicsList::RegisterPhysicsConstructor: <" << name << ">"
           << " fails - name is not defined" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterChemistryConstructor(const G4String& name)
{
  if (name == fChemDNAName)
  {
    return;
  }
  if (verboseLevel > 0)
  {
    G4cout << "===== Register Chemistry constructor ==== " << name << G4endl;
  }
  if (name == "G4EmDNAChemistry")
  {
    delete fEmDNAChemistryList;
    fEmDNAChemistryList = new G4EmDNAChemistry();
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    fChemDNAName = name;
  }
  else if (name == "G4EmDNAChemistry_option1")
  {
    delete fEmDNAChemistryList;
    fEmDNAChemistryList = new G4EmDNAChemistry_option1();
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    fChemDNAName = name;
  }
  else if (name == "G4EmDNAChemistry_option2")
  {
    delete fEmDNAChemistryList;
    fEmDNAChemistryList = new G4EmDNAChemistry_option2();
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    fChemDNAName = name;
  }
  else if (name == "G4EmDNAChemistry_option3")
  {
    delete fEmDNAChemistryList;
    fEmDNAChemistryList = new G4EmDNAChemistry_option3();
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    fChemDNAName = name;
  }
  else if (name == "G4EmDNAChemistry_userOption")
  {
    delete fEmDNAChemistryList;
    fEmDNAChemistryList = new G4EmDNAChemistryUserOption(fDMSO, fOxygen);
    fEmDNAChemistryList->SetVerboseLevel(verboseLevel);
    fChemDNAName = name;
  }
  else
  {
    G4cout << "PhysicsList::RegisterChemistryConstructor: <" << name << ">"
           << " fails - name is not defined" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpPhysicsUI)
  {
    RegisterPhysicsConstructor(newValue);
  }
  if (command == fpChemistryUI)
  {
    RegisterChemistryConstructor(newValue);
  }
  if (command == fpDMSOUI)
  {
    fDMSO = fpDMSOUI->GetNewDoubleValue(newValue);
  }
  else if (command == fpOxygenUI)
  {
    fOxygen = fpOxygenUI->GetNewDoubleValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
