#include "Strand.hh"
#include "Nucleotide.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Nucleotide& Strand::GetNucleotide(G4int index) const {
    // get a nucleotide object from a map
    for (auto& pair : fNucleotides) {
        G4int id = pair.first;
        const Nucleotide& nucleotide = pair.second;
        if (id == index) return pair.second;
    }
    G4cerr << "couldn't find a nucleotide: " << index << G4endl;
    std::exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......