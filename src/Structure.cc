#include "DNAStructure.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid DNAStructure::GetEllipsoid(G4String chainId, G4int residueId, G4String compoundName)
{
  if (fEllipsoidMap.find(DNAComponentKey{chainId, residueId, compoundName}) == fEllipsoidMap.end())
  {
    return Ellipsoid();  // return empty ellipsoid if not found
  }
  // get ellipsoid from the map
  return fEllipsoidMap.at(DNAComponentKey{chainId, residueId, compoundName});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<Plane> DNAStructure::GetPlanes(G4String chainId, G4int residueId,
                                        G4String compoundName)
{
  if (fPlanesMap.find(DNAComponentKey{chainId, residueId, compoundName}) == fPlanesMap.end())
  {
    return std::vector<Plane>();  // return empty vector if not found
  }
  // get planes from the map
  return fPlanesMap.at(DNAComponentKey{chainId, residueId, compoundName});
}

void DNAStructure::SetEllipsoid(G4String chainId, G4int residueId, G4String compoundName,
                             Ellipsoid ellipsoid)
{
  fEllipsoidMap[DNAComponentKey{chainId, residueId, compoundName}] = ellipsoid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::AddPlane(G4String chainId, G4int residueId, G4String compoundName, Plane plane)
{
  fPlanesMap[DNAComponentKey{chainId, residueId, compoundName}].push_back(plane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAStructure::SetPlanes(G4String chainId, G4int residueId, G4String compoundName,
                          std::vector<Plane> planes)
{
  fPlanesMap[DNAComponentKey{chainId, residueId, compoundName}] = planes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......