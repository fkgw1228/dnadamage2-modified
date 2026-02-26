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

#ifndef MODEL_HH
#define MODEL_HH

#include "globals.hh"

#include "Chain.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Model
{
  public:
    // Constructor
    Model() : fId(0) {}
    Model(G4int modelId) : fId(modelId) {}
    Model(G4int modelId, std::map<G4String, Chain> chains) : fId(modelId), fChains(chains) {}

    // Getter
    G4int GetId() { return this->fId; }
    std::map<G4String, Chain>& GetChains() { return this->fChains; }
    Chain& GetChain(G4String chainId);

    // Setter
    void SetId(G4int modelId) { this->fId = modelId; }
    void SetChains(std::map<G4String, Chain> chains) { fChains = chains; }
    void AddChain(Chain chain) { fChains[chain.GetId()] = chain; }

  private:
    G4int fId;  // Model ID
    std::map<G4String, Chain> fChains;  // Chains in the Model keyed by chain ID
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif