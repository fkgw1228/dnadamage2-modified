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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4AffineTransform.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RandomTools.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4VoxelLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fpParticleGun(0)
{
  G4int n_particle = 1;
  fpParticleGun = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  fpParticleGun->SetParticleDefinition(particle);
  fpParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fpParticleGun->SetParticleEnergy(100 * keV);
  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  fpParticleNumUI = new G4UIcmdWithAnInteger("/fpGun/setPrimariesPerEvent", this);
  fpParticleNumUI->SetGuidance("Number of primary particles per event.");

  fpParticleEnergyUI = new G4UIcmdWithADoubleAndUnit("/fpGun/setEnergy", this);
  fpParticleEnergyUI->SetGuidance("Initial energy of primary particles");

  fpIrradiationTypeUI = new G4UIcmdWithAString("/fpGun/setIrradiationType", this);
  fpIrradiationTypeUI->SetGuidance("Set irradiation type (fixed or spherical)");

  fpSphereRadiusUI = new G4UIcmdWithADoubleAndUnit("/fpGun/setSphereRadius", this);
  fpSphereRadiusUI->SetGuidance(
    "Set radius of spherical irradiation volume (used only if irradiation type is spherical)");

  fpParticlePosUI = new G4UIcmdWith3VectorAndUnit("/fpGun/setPosition", this);
  fpParticlePosUI->SetGuidance(
    "Initial position of primary particles (used only if irradiation type is fixed)");

  fpParticleDirUI = new G4UIcmdWith3Vector("/fpGun/setDirection", this);
  fpParticleDirUI->SetGuidance(
    "Initial direction of primary particles (used only if irradiation type is fixed)");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpParticleNumUI;
  delete fpParticleEnergyUI;
  delete fpIrradiationTypeUI;
  delete fpSphereRadiusUI;
  delete fpParticlePosUI;
  delete fpParticleDirUI;
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  for (G4int i = 0; i < fParticleNum; i++)
  {
    if (fIrradiationType == "spherical")
    {
      // random position on sphere
      G4ThreeVector dir = G4RandomDirection();
      fParticlePosition = dir * fSphereRadius;
      // random direction toward center
      fParticleDirection = (-dir + G4RandomDirection() * 0.1).unit();
    }
    fpParticleGun->SetParticleMomentumDirection(fParticleDirection);
    fpParticleGun->SetParticlePosition(fParticlePosition);
    fpParticleGun->SetParticleEnergy(fParticleEnergy);
    fpParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpParticleNumUI)
  {
    fParticleNum = fpParticleNumUI->GetNewIntValue(newValue);
  }
  if (command == fpParticleEnergyUI)
  {
    fParticleEnergy = fpParticleEnergyUI->GetNewDoubleValue(newValue);
  }
  if (command == fpIrradiationTypeUI)
  {
    fIrradiationType = newValue;
  }
  if (command == fpSphereRadiusUI)
  {
    fSphereRadius = fpSphereRadiusUI->GetNewDoubleValue(newValue);
  }
  if (command == fpParticlePosUI)
  {
    fParticlePosition = fpParticlePosUI->GetNew3VectorValue(newValue);
  }
  if (command == fpParticleDirUI)
  {
    fParticleDirection = fpParticleDirUI->GetNew3VectorValue(newValue).unit();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....