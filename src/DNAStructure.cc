#include "DNAStructure.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid DNAStructure::GetEllipsoid(G4String chainId, G4int nucleotideId, G4String compoundName) const {
    if (fEllipsoidMap.find(DNALocationKey{chainId, nucleotideId, compoundName}) == fEllipsoidMap.end()) {
        return Ellipsoid();   // return empty ellipsoid if not found
    }
    // get ellipsoid from the map
    return fEllipsoidMap.at(DNALocationKey{chainId, nucleotideId, compoundName});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<Plane> DNAStructure::GetPlanes(G4String chainId, G4int nucleotideId, G4String compoundName) const {
    if (fPlanesMap.find(DNALocationKey{chainId, nucleotideId, compoundName}) == fPlanesMap.end()) {
        return std::vector<Plane>();   // return empty vector if not found
    }
    // get planes from the map
    return fPlanesMap.at(DNALocationKey{chainId, nucleotideId, compoundName});
}

void DNAStructure::SetEllipsoid(G4String chainId, G4int nucleotideId, G4String compoundName, Ellipsoid ellipsoid) {
    fEllipsoidMap[DNALocationKey{chainId, nucleotideId, compoundName}] = ellipsoid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::AddPlane(G4String chainId, G4int nucleotideId, G4String compoundName, Plane plane) {
    fPlanesMap[DNALocationKey{chainId, nucleotideId, compoundName}].push_back(plane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::SetPlanes(G4String chainId, G4int nucleotideId, G4String compoundName, std::vector<Plane> planes) {
    fPlanesMap[DNALocationKey{chainId, nucleotideId, compoundName}] = planes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......