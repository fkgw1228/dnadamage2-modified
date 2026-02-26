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

#include "DNAMolecules.hh"

#include "G4DNABrownianTransportation.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H1

G4DNADeoxyriboseH1* G4DNADeoxyriboseH1::theInstance = nullptr;
G4DNADeoxyriboseH1* G4DNADeoxyriboseH1::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H1";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H1^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH1*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H2a

G4DNADeoxyriboseH2a* G4DNADeoxyriboseH2a::theInstance = nullptr;
G4DNADeoxyriboseH2a* G4DNADeoxyriboseH2a::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H2a";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H2a^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH2a*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H2b

G4DNADeoxyriboseH2b* G4DNADeoxyriboseH2b::theInstance = nullptr;
G4DNADeoxyriboseH2b* G4DNADeoxyriboseH2b::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H2b";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H2b^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH2b*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H3

G4DNADeoxyriboseH3* G4DNADeoxyriboseH3::theInstance = nullptr;
G4DNADeoxyriboseH3* G4DNADeoxyriboseH3::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H3";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H3^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH3*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H4

G4DNADeoxyriboseH4* G4DNADeoxyriboseH4::theInstance = nullptr;
G4DNADeoxyriboseH4* G4DNADeoxyriboseH4::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H4";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H4^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH4*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H5a

G4DNADeoxyriboseH5a* G4DNADeoxyriboseH5a::theInstance = nullptr;
G4DNADeoxyriboseH5a* G4DNADeoxyriboseH5a::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H5a";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H5a^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH5a*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// Deoxyribose_H5b

G4DNADeoxyriboseH5b* G4DNADeoxyriboseH5b::theInstance = nullptr;
G4DNADeoxyriboseH5b* G4DNADeoxyriboseH5b::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "DNA_Deoxyribose_H5b";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr)
  {
    G4double mass = 1.0079 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0,  // charge
                                          1,  // number of occupancies
                                          0.5 * angstrom, 2);  // radius has to be checked

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0, 1);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName("H5b^0");
  }
  theInstance = static_cast<G4DNADeoxyriboseH5b*>(anInstance);
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Histone

G4DNAHistone* G4DNAHistone::fHistoneInstance = 0;
G4DNAHistone* G4DNAHistone::Definition()
{
  if (fHistoneInstance != 0) return fHistoneInstance;
  const G4String name = "DNA_Histone";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "Histone^{0}";
    G4double mass = 1.0 * g / Avogadro * c_squared;
    anInstance = new G4MoleculeDefinition(name, mass, 1e-150, 0, 0, 1.0 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }

  fHistoneInstance = static_cast<G4DNAHistone*>(anInstance);
  return fHistoneInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H1

G4DNADamagedDeoxyriboseH1OH* G4DNADamagedDeoxyriboseH1OH::fDamagedDeoxyriboseH1OHInstance = 0;
G4DNADamagedDeoxyriboseH1OH* G4DNADamagedDeoxyriboseH1OH::Definition()
{
  if (fDamagedDeoxyriboseH1OHInstance != 0) return fDamagedDeoxyriboseH1OHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H1_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H1_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH1OHInstance = static_cast<G4DNADamagedDeoxyriboseH1OH*>(anInstance);
  return fDamagedDeoxyriboseH1OHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H2a

G4DNADamagedDeoxyriboseH2aOH* G4DNADamagedDeoxyriboseH2aOH::fDamagedDeoxyriboseH2aOHInstance = 0;
G4DNADamagedDeoxyriboseH2aOH* G4DNADamagedDeoxyriboseH2aOH::Definition()
{
  if (fDamagedDeoxyriboseH2aOHInstance != 0) return fDamagedDeoxyriboseH2aOHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H2a_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H2a_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH2aOHInstance = static_cast<G4DNADamagedDeoxyriboseH2aOH*>(anInstance);
  return fDamagedDeoxyriboseH2aOHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H2b

G4DNADamagedDeoxyriboseH2bOH* G4DNADamagedDeoxyriboseH2bOH::fDamagedDeoxyriboseH2bOHInstance = 0;
G4DNADamagedDeoxyriboseH2bOH* G4DNADamagedDeoxyriboseH2bOH::Definition()
{
  if (fDamagedDeoxyriboseH2bOHInstance != 0) return fDamagedDeoxyriboseH2bOHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H2b_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H2b_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH2bOHInstance = static_cast<G4DNADamagedDeoxyriboseH2bOH*>(anInstance);
  return fDamagedDeoxyriboseH2bOHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H3

G4DNADamagedDeoxyriboseH3OH* G4DNADamagedDeoxyriboseH3OH::fDamagedDeoxyriboseH3OHInstance = 0;
G4DNADamagedDeoxyriboseH3OH* G4DNADamagedDeoxyriboseH3OH::Definition()
{
  if (fDamagedDeoxyriboseH3OHInstance != 0) return fDamagedDeoxyriboseH3OHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H3_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H3_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH3OHInstance = static_cast<G4DNADamagedDeoxyriboseH3OH*>(anInstance);
  return fDamagedDeoxyriboseH3OHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H4

G4DNADamagedDeoxyriboseH4OH* G4DNADamagedDeoxyriboseH4OH::fDamagedDeoxyriboseH4OHInstance = 0;
G4DNADamagedDeoxyriboseH4OH* G4DNADamagedDeoxyriboseH4OH::Definition()
{
  if (fDamagedDeoxyriboseH4OHInstance != 0) return fDamagedDeoxyriboseH4OHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H4_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H4_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH4OHInstance = static_cast<G4DNADamagedDeoxyriboseH4OH*>(anInstance);
  return fDamagedDeoxyriboseH4OHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H5a

G4DNADamagedDeoxyriboseH5aOH* G4DNADamagedDeoxyriboseH5aOH::fDamagedDeoxyriboseH5aOHInstance = 0;
G4DNADamagedDeoxyriboseH5aOH* G4DNADamagedDeoxyriboseH5aOH::Definition()
{
  if (fDamagedDeoxyriboseH5aOHInstance != 0) return fDamagedDeoxyriboseH5aOHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H5a_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H5a_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH5aOHInstance = static_cast<G4DNADamagedDeoxyriboseH5aOH*>(anInstance);
  return fDamagedDeoxyriboseH5aOHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
// OH_Damamged_Deoxyribose_H5b

G4DNADamagedDeoxyriboseH5bOH* G4DNADamagedDeoxyriboseH5bOH::fDamagedDeoxyriboseH5bOHInstance = 0;
G4DNADamagedDeoxyriboseH5bOH* G4DNADamagedDeoxyriboseH5bOH::Definition()
{
  if (fDamagedDeoxyriboseH5bOHInstance != 0) return fDamagedDeoxyriboseH5bOHInstance;
  const G4String name = "DNA_DamagedDeoxyribose_H5b_OH";

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == 0)
  {
    const G4String formatedName = "DamagedDeoxyribose_H5b_OH^{0}";

    G4double mass = 31.99546 * g / Avogadro * c_squared;
    anInstance =
      new G4MoleculeDefinition(name, mass, 1e-150 * (m * m / s), 0, 0, 1.7 * angstrom, 2);

    ((G4MoleculeDefinition*)anInstance)->SetLevelOccupation(0);
    ((G4MoleculeDefinition*)anInstance)->SetFormatedName(formatedName);
  }
  fDamagedDeoxyriboseH5bOHInstance = static_cast<G4DNADamagedDeoxyriboseH5bOH*>(anInstance);
  return fDamagedDeoxyriboseH5bOHInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....