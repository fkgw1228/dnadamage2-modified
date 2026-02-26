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

#ifndef RESIDUE_HH
#define RESIDUE_HH

#include "globals.hh"

#include "Compound.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Residue
{
  public:
    // Constructor
    Residue() = default;
    //Residue(G4int residueId) : fId(residueId), fName("") {}
    Residue(G4int residueId, G4String name) : fId(residueId), fName(name) {}

    // Getter
    G4int GetId() { return this->fId; }
    G4String GetName() { return this->fName; }
    std::map<G4String, Compound>& GetCompounds() { return this->fCompounds; }
    Compound& GetCompound(G4String compoundName);

    // Setter
    void SetName(G4String name) { this->fName = name; }
    void AddCompound(Compound compound) { this->fCompounds[compound.GetName()] = compound; }

  private:
    G4int fId;                                // residue number
    G4String fName;                           // residue name (e.g., DA,DC,DC,DT,)
    std::map<G4String, Compound> fCompounds;  // pairs of name and compound constituting the residue
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif