#ifndef MODEL_HH
#define MODEL_HH

#include "globals.hh"

#include "Chain.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Model
{
  public:
    // Constructor
    Model() : fId(0) {}
    Model(G4int modelId) : fId(modelId) {}
    Model(G4int modelId, std::map<G4String, Chain> chains) : fId(modelId), fChains(chains) {}

    // Getter
    G4int GetId() { return this->fId; }
    std::map<G4String, Chain>& GetChains() { return this->fChains; }
    Chain& GetChain(G4String chainId);

    // Setter
    void SetId(G4int modelId) { this->fId = modelId; }
    void SetChains(std::map<G4String, Chain> chains) { fChains = chains; }
    void AddChain(Chain chain) { fChains[chain.GetId()] = chain; }

  private:
    G4int fId;  // Model ID
    std::map<G4String, Chain> fChains;  // Chains in the Model keyed by chain ID
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif