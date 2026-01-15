#ifndef DNA_HH
#define DNA_HH

#include <map>
#include "globals.hh"
#include "Strand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNA {
public:
    // Constructor
    DNA(): fId(0) {}
    DNA(G4int dnaId): fId(dnaId) {}
    DNA(G4int dnaId, std::map<G4String, Strand> strands)
        : fId(dnaId), fStrands(strands) {}

    // Getter
    G4int GetId() const { return this->fId; }
    const std::map<G4String, Strand>& GetStrands() const { return this->fStrands; }
    const Strand& GetStrand(G4String chainId) const;

    // Setter
    void SetId(G4int dnaId) { this->fId = dnaId; }
    void SetStrands(std::map<G4String, Strand> strands) { fStrands = strands; }
    void AddStrand(Strand strand) { fStrands[strand.GetId()] = strand; }

private:
    G4int fId;                               // DNA ID
    std::map<G4String, Strand> fStrands;     // Strands in the DNA keyed by chain ID
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif