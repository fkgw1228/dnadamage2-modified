#ifndef DNASTRUCTURE_HH
#define DNASTRUCTURE_HH

#include "Ellipsoid.hh"
#include "Model.hh"
#include "Plane.hh"

#include <map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct DNAComponentKey
{
    G4String chainId;
    G4int residueId;
    G4String compoundName;

    // Define less-than operator for use in std::map
    bool operator<(const DNAComponentKey& other) const
    {
      if (chainId != other.chainId) return chainId < other.chainId;
      if (residueId != other.residueId) return residueId < other.residueId;
      return compoundName < other.compoundName;
    }
};

typedef std::map<DNAComponentKey, Ellipsoid> DNAEllipsoidMap;
typedef std::map<DNAComponentKey, std::vector<Plane>> DNAPlanesMap;

class DNAStructure
{
  public:
    // Constructor
    DNAStructure() {}
    DNAStructure(Model model) : fModel(model) {}

    // Getter
    Model& GetModel() { return fModel; }
    Ellipsoid GetEllipsoid(G4String chainId, G4int residueId, G4String compoundName);
    std::vector<Plane> GetPlanes(G4String chainId, G4int residueId, G4String compoundName);

    // Setter
    void SetModel(Model model) { fModel = model; }
    void SetEllipsoid(G4String chainId, G4int residueId, G4String compoundName, Ellipsoid ellipsoid);
    void AddPlane(G4String chainId, G4int residueId, G4String compoundName, Plane plane);
    void SetPlanes(G4String chainId, G4int residueId, G4String compoundName,
                   std::vector<Plane> planes);

  private:
    Model fModel;
    DNAEllipsoidMap fEllipsoidMap;
    DNAPlanesMap fPlanesMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif