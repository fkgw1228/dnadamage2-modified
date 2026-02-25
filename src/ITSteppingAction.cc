#include "ITSteppingAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITSteppingAction::ITSteppingAction() : G4UserSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITSteppingAction::~ITSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* track = aStep->GetTrack();
  const G4String& moleculeName = track->GetDefinition()->GetParticleName();

  G4bool isDNAMolecule = G4StrUtil::contains(moleculeName, "DNA_Deoxyribose")
                         || G4StrUtil::contains(moleculeName, "DNA_Histone");
  if (isDNAMolecule)
  {
    return;
  }

  const G4String& volumeName = preStepPoint->GetPhysicalVolume()->GetName();
  G4bool isInHistone = false;
  for (const G4String& aminoAcidName : fAminoAcidNames)
  {
    if (G4StrUtil::contains(volumeName, aminoAcidName))
    {
      isInHistone = true;
      break;
    }
  }
  // remove non-DNA molecules inside DNA and histone volumes
  if (isInHistone && !isDNAMolecule)
  {
    track->SetTrackStatus(fStopAndKill);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
