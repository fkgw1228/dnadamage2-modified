#include "Chain.hh"

#include "Residue.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue& Chain::GetResidue(G4int residueId)
{
  // get a residue object from a map
  for (auto& [id, residue] : fResidues)
  {
    if (id == residueId) return residue;
  }
  G4cerr << "couldn't find a residue: " << residueId << G4endl;
  std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......