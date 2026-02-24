#include "Residue.hh"

#include "Compound.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Compound& Residue::GetCompound(G4String compoundName)
{
  // get a compound object from a map
  for (auto& [name, compound] : fCompounds)
  {
    if (name == compoundName) return compound;
  }
  G4cerr << "couldn't find a compound: " << compoundName << G4endl;
  std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......