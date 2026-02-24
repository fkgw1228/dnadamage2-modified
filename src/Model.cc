#include "Model.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Chain& Model::GetChain(G4String chainId)
{
  // get a chain object from a map
  for (auto& [id, chain] : fChains)
  {
    if (id == chainId) return chain;
  }
  G4cerr << "couldn't find a chain: " << chainId << G4endl;
  std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......