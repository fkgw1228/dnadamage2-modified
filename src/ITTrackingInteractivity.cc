//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Original license: Geant4 Software License (see GEANT4_LICENSE).
// Copyright (C) 2026 Shun Fukagawa, Tsukasa Aso

#include "ITTrackingInteractivity.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4IT.hh"
#include "G4RichTrajectory.hh"
#include "G4SmoothTrajectory.hh"
#include "G4TrackingInformation.hh"
#include "G4Trajectory.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VSteppingVerbose.hh"
#include "G4VTrajectory.hh"
#include "G4VisManager.hh"

ITTrackingInteractivity::ITTrackingInteractivity()
{
  fpUserTrackingAction = 0;
  fpUserSteppingAction = 0;

  ////////////////////////////
  // In case you want to use same tracking/stepping action
  // for normal and IT stepping
  /*
  fpUserTrackingAction =
  trackingManager->GetUserTrackingAction();
  fpUserSteppingAction =
  G4EventManager::GetEventManager()->GetUserSteppingAction();
 */
  ////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITTrackingInteractivity::~ITTrackingInteractivity()
{
  G4EventManager* eventManager = G4EventManager::GetEventManager();

  if (eventManager)
  {
    G4UserTrackingAction* std_trackAct = eventManager->GetUserTrackingAction();
    if (fpUserTrackingAction != std_trackAct && fpUserTrackingAction) delete fpUserTrackingAction;

    G4UserSteppingAction* std_stepAct = eventManager->GetUserSteppingAction();
    if (fpUserSteppingAction != std_stepAct && fpUserSteppingAction) delete fpUserSteppingAction;
  }
  else
  {
    if (fpUserSteppingAction)
    {
      delete fpUserSteppingAction;
    }
    if (fpUserTrackingAction)
    {
      delete fpUserTrackingAction;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::Initialize()
{
  G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::StartTracking(G4Track* track)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::AppendStep(G4Track* track, G4Step* step)
{
  ////////////////////////////
  // If you want to use sensitive detector
  /*
    G4VPhysicalVolume* currentVolume =
    step->GetPreStepPoint()->GetPhysicalVolume();
    G4SteppingControl stepControlFlag =  step->GetControlFlag();

    if( currentVolume  != 0 && stepControlFlag != AvoidHitInvocation)
    {
        G4VSensitiveDetector* sensitive = step->GetPreStepPoint()->
                GetSensitiveDetector();
        if( sensitive != 0 )
        {
            sensitive->Hit(fpStep);
        }
    }
   */
  ////////////////////////////

  if (fpUserSteppingAction) fpUserSteppingAction->UserSteppingAction(step);

  ////////////////////////////
  // If you want to use regional stepping action
  /*
  G4UserSteppingAction* regionalAction
          = fpStep->GetPreStepPoint()->GetPhysicalVolume()->
                GetLogicalVolume()->GetRegion()->
                GetRegionalSteppingAction();
  if( regionalAction ) regionalAction->UserSteppingAction(fpStep);
  */
  ////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::EndTracking(G4Track* track)
{
  //     // Post tracking user intervention process.
  //     if (fpUserTrackingAction != 0) {
  //         fpUserTrackingAction->PostUserTrackingAction(track);
  //     }
}

void ITTrackingInteractivity::Finalize()
{
}
