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

#ifndef ATOM_HH
#define ATOM_HH

#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Atom
{
  public:
    // Constructor
    Atom(G4String name, G4int id, G4String element, G4ThreeVector position, G4double radius)
      : fName(name), fId(id), fElement(element), fPosition(position), fRadius(radius) {}

    // Getters
    G4String GetName() { return fName; }
    G4int GetId() { return fId; }
    G4String GetElement() { return fElement; }
    G4ThreeVector GetPosition() { return fPosition; }
    G4double GetRadius() { return fRadius; }

    // Setter
    void SetPosition(G4ThreeVector newPosition) { fPosition = newPosition; }

  private:
    G4String fName;           // Atom name
    G4int fId;                // Atom ID
    G4String fElement;        // Element
    G4ThreeVector fPosition;  // Atomic position
    G4double fRadius;         // Atom's radius
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif