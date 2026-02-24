#ifndef DNAHANDLER_HH
#define DNAHANDLER_HH

#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include "Model.hh"
#include "Structure.hh"

#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// * To use this class, call its static methods like DNAHandler::GetMovedDNA(...)
class DNAHandler
{
  public:
    // 　Move the entire DNA structure by a given offset
    static Structure GetMovedStructure(G4int, Structure&, G4ThreeVector);

    // 　Get the center coordinates of the DNA
    static G4ThreeVector GetCenter(Model& model);

    // 　Get the range (minimum and maximum coordinates) of the DNA
    static std::pair<G4ThreeVector, G4ThreeVector> GetRange(Model& model);

    // 　Print the information of the DNA structure to the console
    static void Print(Structure& structure);

  private:
    DNAHandler() = delete;  // Prevent instantiation of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif