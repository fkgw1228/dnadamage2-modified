
#ifndef ITTRACKINGINTERACTIVITY_HH
#define ITTRACKINGINTERACTIVITY_HH 1

#include "G4ITTrackingInteractivity.hh"

#include <vector>

class G4VTrajectory;

class ITTrackingInteractivity : public G4ITTrackingInteractivity
{
    G4UserTrackingAction* fpUserTrackingAction;
    G4UserSteppingAction* fpUserSteppingAction;
    // int fStoreTrajectory;
    // std::vector<G4VTrajectory*> fTrajectories;

  public:
    ITTrackingInteractivity();
    virtual ~ITTrackingInteractivity();

    virtual void Initialize();
    virtual void StartTracking(G4Track*);
    virtual void AppendStep(G4Track* track, G4Step* step);
    virtual void EndTracking(G4Track*);
    virtual void Finalize();

    void SetUserAction(G4UserTrackingAction* trackAct) { fpUserTrackingAction = trackAct; }
    inline G4UserTrackingAction* GetUserTrackingAction() { return fpUserTrackingAction; }

    void SetUserAction(G4UserSteppingAction* stepAct) { fpUserSteppingAction = stepAct; }
    inline G4UserSteppingAction* GetUserSteppingAction() { return fpUserSteppingAction; }
};

#endif