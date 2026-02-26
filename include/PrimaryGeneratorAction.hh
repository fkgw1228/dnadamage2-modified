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

#ifndef DNADAMAGE2_PrimaryGeneratorAction_h
#define DNADAMAGE2_PrimaryGeneratorAction_h 1

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <G4UImessenger.hh>

#include <cmath>

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event*) override;
    const G4ParticleGun* GetParticleGun() const { return fpParticleGun; }
    void SetNewValue(G4UIcommand* command, G4String newValue) override;

  private:
    // command-line UI
    G4ParticleGun* fpParticleGun = nullptr;
    G4UIcmdWithADoubleAndUnit* fpParticleEnergyUI = nullptr;
    G4UIcmdWithAnInteger* fpParticleNumUI = nullptr;

    G4UIcmdWithAString* fpIrradiationTypeUI = nullptr;
    G4UIcmdWithADoubleAndUnit* fpSphereRadiusUI = nullptr;

    G4UIcmdWith3VectorAndUnit* fpParticlePosUI = nullptr;
    G4UIcmdWith3Vector* fpParticleDirUI = nullptr;

    G4int fParticleNum = 1;
    G4double fParticleEnergy = 100 * keV;
    G4String fIrradiationType = "fixed";
    G4double fSphereRadius = 20 * nanometer;
    G4ThreeVector fParticlePosition = G4ThreeVector(0, 0, 0) * angstrom;
    G4ThreeVector fParticleDirection = G4ThreeVector(1, 0, 0);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
