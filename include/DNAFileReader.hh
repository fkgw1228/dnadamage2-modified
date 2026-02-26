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

#ifndef DNAREADER_HH
#define DNAREADER_HH

#include "G4ThreeVector.hh"

#include "Atom.hh"
#include "Chain.hh"
#include "Compound.hh"
#include "Ellipsoid.hh"
#include "Model.hh"
#include "Plane.hh"
#include "Residue.hh"
#include "DNAStructure.hh"
#include <nlohmann/json.hpp>

#include <string>
#include <vector>

using json = nlohmann::json;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAFileReader
{
  public:
    // Constructor
    DNAFileReader(G4String fileName) : fFileName(fileName) {}

    // Setter
    void SetFileName(G4String fileName) { this->fFileName = fileName; }

    DNAStructure ReadDNAFile(G4int modelId);

  private:
    Model GenerateModel(G4int modelId, const json& dnaInfo);
    Chain GenerateChain(G4String chainId, const json& chainInfo);
    Residue GenerateResidue(G4int residueId, const json& residueInfo);
    Compound GenerateCompound(G4String compoundName, const json& compoundInfo);
    Atom GenerateAtom(const json& atomInfo);
    DNAStructure GenerateStructure(const json& dnaInfo);
    Ellipsoid GenerateEllipsoid(const json& ellipsoidInfo);
    Plane GeneratePlane(const json& planeInfo);

    G4String fFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif