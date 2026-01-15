#ifndef NUCLEOTIDE_HH
#define NUCLEOTIDE_HH

#include  <map>
#include "globals.hh"
#include "Compound.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Nucleotide {
public:
    // Constructor
    Nucleotide() : fId(0), fBaseType("") {}
    Nucleotide(G4int nucleotideId) : fId(nucleotideId), fBaseType("") {}
    Nucleotide(G4int nucleotideId, G4String baseType)
               : fId(nucleotideId), fBaseType(baseType) {}

    // Getter
    G4int GetId() const { return this->fId; }
    G4String GetBaseType() const { return this->fBaseType; }
    const std::map<G4String, Compound>& GetCompounds() const { return this->fCompounds; }
    const Compound& GetCompound(G4String compoundName) const;

    // Setter
    void SetBaseType(G4String baseType) { this->fBaseType = baseType; }
    void AddCompound(Compound compound) { this->fCompounds[compound.GetName()] = compound; }
private:
    G4int fId;                            // nucleotide number
    G4String fBaseType;                   // base_type (e.g., DA,DC,DC,DT)
    std::map<G4String, Compound> fCompounds;   // pairs of name and compound constituting the nucleotide
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif