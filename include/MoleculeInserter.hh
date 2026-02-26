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

#ifndef DNADAMAGE2_moleculeinserter_h
#define DNADAMAGE2_moleculeinserter_h 1

#include "G4ITGun.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <G4memory.hh>

#include <map>
#include <vector>

class G4Track;
class MoleculeInserter;
class G4Molecule;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MoleculeShoot
{
  public:
    MoleculeShoot();
    ~MoleculeShoot() = default;
    void Shoot(MoleculeInserter*);

  public:
    G4int fNumber = 0;
    G4String fMoleculeName = "None";
    G4double fTime = 1;
    G4ThreeVector fPosition = G4ThreeVector();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MoleculeInserter : public G4ITGun
{
  public:
    MoleculeInserter();
    MoleculeInserter(G4bool);
    ~MoleculeInserter() = default;
    void DefineTracks() override;
    void AddMolecule(const G4String& moleculeName, const G4ThreeVector& position, double time = 0);
    void PushToChemistry(G4Molecule*, G4double, G4ThreeVector);
    void CreateMolecule(G4Molecule*, G4double, G4ThreeVector);
    void CreateMolecule(G4String, G4double, G4ThreeVector);
    void Clean();

    std::vector<G4Track*> const GetInsertedTracks() { return fInsertedTracks; }

  private:
    G4bool fSaveTrackID = true;
    std::vector<G4Track*> fInsertedTracks;

  protected:
    std::vector<MoleculeShoot> fShoots;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
