#ifndef STRAND_HH
#define STRAND_HH

#include <vector>
#include <map>
#include "globals.hh"
#include "Nucleotide.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Strand {
public:
    // Constructor
    Strand() : fId("") {}
    Strand(G4String id) : fId(id) {}

    // Getter
    G4String GetId() const { return fId; }
    const std::map<G4int, Nucleotide>& GetNucleotides() const { return fNucleotides; }
    const Nucleotide& GetNucleotide(G4int nucleotideId) const;

    // Setter
    void AddNucleotide(Nucleotide nucleotide) { fNucleotides[nucleotide.GetId()] = nucleotide; }
    
private:
    G4String fId;
    std::map<G4int, Nucleotide> fNucleotides;
};

#endif