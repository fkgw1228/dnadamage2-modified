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
// * tecH1ical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// This example is provided by the Geant4-DNA collaboration
// DNADAMAGE2 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//      J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 
//
/// \file PlasmidMolecules.hh
/// \brief Definition of the additional Plasmid DNA molecules

#ifndef DNADAMAGE2_DNAMolecules_h
#define DNADAMAGE2_DNAMolecules_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MoleculeDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH1 : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH1* theInstance;
    G4DNADeoxyriboseH1() {}
    ~G4DNADeoxyriboseH1() override = default;

public:
    static G4DNADeoxyriboseH1* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH2a : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH2a* theInstance;
    G4DNADeoxyriboseH2a() {}
    ~G4DNADeoxyriboseH2a() override = default;

public:
    static G4DNADeoxyriboseH2a* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH2b : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH2b* theInstance;
    G4DNADeoxyriboseH2b() {}
    ~G4DNADeoxyriboseH2b() override = default;

public:
    static G4DNADeoxyriboseH2b* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH3 : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH3* theInstance;
    G4DNADeoxyriboseH3() {}
    ~G4DNADeoxyriboseH3() override = default;

public:
    static G4DNADeoxyriboseH3* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH4 : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH4* theInstance;
    G4DNADeoxyriboseH4() {}
    ~G4DNADeoxyriboseH4() override = default;

public:
    static G4DNADeoxyriboseH4* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH5a : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH5a* theInstance;
    G4DNADeoxyriboseH5a() {}
    ~G4DNADeoxyriboseH5a() override = default;

public:
    static G4DNADeoxyriboseH5a* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyriboseH5b : public G4MoleculeDefinition
{
private:
    static /*G4ThreadLocal*/ G4DNADeoxyriboseH5b* theInstance;
    G4DNADeoxyriboseH5b() {}
    ~G4DNADeoxyriboseH5b() override = default;

public:
    static G4DNADeoxyriboseH5b* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADeoxyribose : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADeoxyribose* fDeoxyriboseInstance;
  G4DNADeoxyribose() {}
  ~G4DNADeoxyribose() override  = default;

public:
  static G4DNADeoxyribose* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseOH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseOH* fDamagedDeoxyriboseOHInstance;
  G4DNADamagedDeoxyriboseOH() {}
  ~G4DNADamagedDeoxyriboseOH() override = default;

public:
  static G4DNADamagedDeoxyriboseOH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH1OH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH1OH* fDamagedDeoxyriboseH1OHInstance;
  G4DNADamagedDeoxyriboseH1OH() {}
  ~G4DNADamagedDeoxyriboseH1OH() override = default;

public:
  static G4DNADamagedDeoxyriboseH1OH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH2aOH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH2aOH* fDamagedDeoxyriboseH2aOHInstance;
  G4DNADamagedDeoxyriboseH2aOH() {}
  ~G4DNADamagedDeoxyriboseH2aOH() override = default;

public:
  static G4DNADamagedDeoxyriboseH2aOH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH2bOH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH2bOH* fDamagedDeoxyriboseH2bOHInstance;
  G4DNADamagedDeoxyriboseH2bOH() {}
  ~G4DNADamagedDeoxyriboseH2bOH() override = default;

public:
  static G4DNADamagedDeoxyriboseH2bOH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH3OH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH3OH* fDamagedDeoxyriboseH3OHInstance;
  G4DNADamagedDeoxyriboseH3OH() {}
  ~G4DNADamagedDeoxyriboseH3OH() override = default;

public:
  static G4DNADamagedDeoxyriboseH3OH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH4OH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH4OH* fDamagedDeoxyriboseH4OHInstance;
  G4DNADamagedDeoxyriboseH4OH() {}
  ~G4DNADamagedDeoxyriboseH4OH() override = default;

public:
  static G4DNADamagedDeoxyriboseH4OH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH5aOH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH5aOH* fDamagedDeoxyriboseH5aOHInstance;
  G4DNADamagedDeoxyriboseH5aOH() {}
  ~G4DNADamagedDeoxyriboseH5aOH() override = default;

public:
  static G4DNADamagedDeoxyriboseH5aOH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DNADamagedDeoxyriboseH5bOH : public G4MoleculeDefinition
{
private:
  static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH5bOH* fDamagedDeoxyriboseH5bOHInstance;
  G4DNADamagedDeoxyriboseH5bOH() {}
  ~G4DNADamagedDeoxyriboseH5bOH() override = default;

public:
  static G4DNADamagedDeoxyriboseH5bOH* Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// class G4DNADamagedDeoxyriboseH : public G4MoleculeDefinition
// {
// private:
//   static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseH* fDamagedDeoxyriboseHInstance;
//   G4DNADamagedDeoxyriboseH() {}
//   ~G4DNADamagedDeoxyriboseH() override  = default;

// public:
//   static G4DNADamagedDeoxyriboseH* Definition();
// };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// class G4DNADamagedDeoxyriboseEAQ : public G4MoleculeDefinition
// {
// private:
//   static /*G4ThreadLocal*/ G4DNADamagedDeoxyriboseEAQ* fDamagedDeoxyriboseEAQInstance;
//   G4DNADamagedDeoxyriboseEAQ() {}
//   ~G4DNADamagedDeoxyriboseEAQ() override  = default;

// public:
//   static G4DNADamagedDeoxyriboseEAQ* Definition();
// };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
