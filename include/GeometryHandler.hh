#ifndef GeometryHandler_HH
#define GeometryHandler_HH

#include "G4RotationMatrix.hh"

#include "Ellipsoid.hh"
#include "Plane.hh"

#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GeometryHandler
{
  public:
    static G4RotationMatrix* GetPlaneRotationMatrix(Plane plane);
    static G4ThreeVector GetPlaneTranslationForEllipsoid(Ellipsoid ellipsoid, Plane plane);

  private:
    GeometryHandler() = delete;  // Prevent instantiation of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif