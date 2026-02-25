#include "G4ChemTimeStepModel.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include "globals.hh"

#ifndef G4EmDNAChemistryUserOption2_hh
#define G4EmDNAChemistryUserOption2_hh 1

class G4DNAMolecularReactionTable;
class G4DNAMolecularReactionData;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmDNAChemistryUserOption : public G4VUserChemistryList, public G4VPhysicsConstructor
{
  public:
    G4EmDNAChemistryUserOption();
    G4EmDNAChemistryUserOption(G4double dmso, G4double oxygen);
    ~G4EmDNAChemistryUserOption() override = default;

    void ConstructParticle() override { ConstructMolecule(); }
    void ConstructMolecule() override;
    void ConstructProcess() override;

    void ConstructDissociationChannels() override;
    void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable) override;
    void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable) override;

  private:
    static void SetReactionType(G4DNAMolecularReactionData* pData, G4ChemTimeStepModel model);

    G4double fDMSO = 0.;
    G4double fOxygen = 0.;
};

#endif

//...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......