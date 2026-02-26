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
//
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

#ifndef DNASTRUCTURE_HH
#define DNASTRUCTURE_HH

#include "Ellipsoid.hh"
#include "Model.hh"
#include "Plane.hh"

#include <map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct DNAComponentKey
{
    G4String chainId;
    G4int residueId;
    G4String compoundName;

    // Define less-than operator for use in std::map
    bool operator<(const DNAComponentKey& other) const
    {
      if (chainId != other.chainId) return chainId < other.chainId;
      if (residueId != other.residueId) return residueId < other.residueId;
      return compoundName < other.compoundName;
    }
};

typedef std::map<DNAComponentKey, Ellipsoid> DNAEllipsoidMap;
typedef std::map<DNAComponentKey, std::vector<Plane>> DNAPlanesMap;

class DNAStructure
{
  public:
    // Constructor
    DNAStructure() {}
    DNAStructure(Model model) : fModel(model) {}

    // Getter
    Model& GetModel() { return fModel; }
    Ellipsoid GetEllipsoid(G4String chainId, G4int residueId, G4String compoundName);
    std::vector<Plane> GetPlanes(G4String chainId, G4int residueId, G4String compoundName);

    // Setter
    void SetModel(Model model) { fModel = model; }
    void SetEllipsoid(G4String chainId, G4int residueId, G4String compoundName, Ellipsoid ellipsoid);
    void AddPlane(G4String chainId, G4int residueId, G4String compoundName, Plane plane);
    void SetPlanes(G4String chainId, G4int residueId, G4String compoundName,
                   std::vector<Plane> planes);

  private:
    Model fModel;
    DNAEllipsoidMap fEllipsoidMap;
    DNAPlanesMap fPlanesMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif