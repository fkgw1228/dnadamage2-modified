#ifndef DNAHANDLER_HH
#define DNAHANDLER_HH

#include <utility>
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "DNAStructure.hh"
#include "DNA.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAHandler
{
public:
    //　Move the entire DNA structure by a given offset
    static DNAStructure GetMovedDNA(G4int dnaId,
                                    const DNAStructure& dnaStructure,
                                    const G4ThreeVector& offset);

    //　Get the center coordinates of the DNA
    static G4ThreeVector GetCenter(const DNA& dna);

    //　Get the range (minimum and maximum coordinates) of the DNA
    static std::pair<G4ThreeVector, G4ThreeVector> GetRange(const DNA& dna);

    //　Print the information of the DNA structure to the console
    static void Print(const DNAStructure& dnaStructure);

private:
    DNAHandler() = delete;  // Prevent instantiation of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif