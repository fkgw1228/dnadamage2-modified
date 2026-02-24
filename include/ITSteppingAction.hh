#ifndef ITSTEPPINGACTION_HH
#define ITSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <vector>

class ITSteppingAction : public G4UserSteppingAction
{
  public:
    ITSteppingAction();
    virtual ~ITSteppingAction();
    virtual void UserSteppingAction(const G4Step*);

  private:
    std::vector<G4String> fAminoAcidNames = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                             "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
};

#endif