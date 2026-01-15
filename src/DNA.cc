#include "DNA.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Strand& DNA::GetStrand(G4String chainId) const {
    // get a strand object from a map
    for (auto& pair : fStrands) {
        G4String id = pair.first;
        Strand strand = pair.second;
        if (id == chainId) return pair.second;
    }
    G4cerr << "couldn't find a strand: " << chainId << G4endl;
    std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......