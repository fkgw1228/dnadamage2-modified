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

#ifndef DNAHANDLER_HH
#define DNAHANDLER_HH

#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include "Model.hh"
#include "DNAStructure.hh"

#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// * To use this class, call its static methods like DNAHandler::GetMovedDNA(...)
class DNAHandler
{
  public:
    // 　Move the entire DNA structure by a given offset
    static DNAStructure GetMovedStructure(G4int, DNAStructure&, G4ThreeVector);

    // 　Get the center coordinates of the DNA
    static G4ThreeVector GetCenter(Model& model);

    // 　Get the range (minimum and maximum coordinates) of the DNA
    static std::pair<G4ThreeVector, G4ThreeVector> GetRange(Model& model);

    // 　Print the information of the DNA structure to the console
    static void Print(DNAStructure& structure);

  private:
    DNAHandler() = delete;  // Prevent instantiation of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif