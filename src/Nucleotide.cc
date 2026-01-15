#include "Nucleotide.hh"
#include "Compound.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Compound& Nucleotide::GetCompound(G4String compoundName) const {
    // get a compound object from a map
    for(const auto& pair : fCompounds) {
        G4String name = pair.first;
        const Compound& compound = pair.second;
        if(name == compoundName) return pair.second;
    }
    G4cerr << "couldn't find a compound: " << compoundName << G4endl;
    std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......