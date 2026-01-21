#include "GeometryHandler.hh"
#include "G4INCLGlobals.hh"
#include "G4SystemOfUnits.hh"
#include "Eigen/Dense"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix* GeometryHandler::GetPlaneRotationMatrix(Plane plane) {
    /* Obtain the rotation matrix that transforms the xy plane at z=0 to the specified plane */
    G4ThreeVector w = plane.GetW().unit();

    // rotate xy plane and conform upper side of box with true plane
    G4ThreeVector n(0.0, 0.0, 1.0);  // normal vector of xy plane
    
    // Handle special cases when n and w are parallel or anti-parallel
    G4double dotProduct = n.dot(w);
    G4RotationMatrix* planeRotation = new G4RotationMatrix();
    
    if (std::abs(dotProduct - 1.0) < 1e-10) {
        // n and w are parallel (same direction) - no rotation needed
        return planeRotation;  // identity matrix
    } else if (std::abs(dotProduct + 1.0) < 1e-10) {
        // n and w are anti-parallel (opposite direction) - rotate 180° around x-axis
        planeRotation->rotateX(CLHEP::pi);
        return planeRotation;
    }
    G4double theta = G4INCL::Math::arcCos(dotProduct);  // angle between the plane and xy plane
    // rotation axis
    G4ThreeVector rotationAxis = (n.cross(w)).unit();   // cross product
    G4double n1 = rotationAxis.x(), n2 = rotationAxis.y(), n3 = rotationAxis.z();
    
    // Rodrigues' rotation formula
    planeRotation->setRows(
        G4ThreeVector(n1*n1*(1-cos(theta))+   cos(theta), n1*n2*(1-cos(theta))+n3*sin(theta), n1*n3*(1-cos(theta))-n2*sin(theta)),
        G4ThreeVector(n2*n1*(1-cos(theta))-n3*sin(theta), n2*n2*(1-cos(theta))+   cos(theta), n2*n3*(1-cos(theta))+n1*sin(theta)),
        G4ThreeVector(n1*n3*(1-cos(theta))+n2*sin(theta), n2*n3*(1-cos(theta))-n1*sin(theta), n3*n3*(1-cos(theta))+   cos(theta))
    );
    return planeRotation;
}

G4ThreeVector GeometryHandler::GetPlaneTranslationForEllipsoid(Ellipsoid ellipsoid, Plane plane) {
    /* Return the translation vector from the origin to the center of the intersection
       between the ellipsoid and the plane, The vector is given by:
       c - (b + wc) / (wE^-1w) * (E^-1 w)
       where c is the ellipsoid center, b is the plane bias, w is the plane normal vector,
       and E is the ellipsoid matrix.
    */
    Eigen::Matrix3d E = ellipsoid.GetEMatrix();
    G4ThreeVector c = ellipsoid.GetCenter();
    G4ThreeVector w = plane.GetW().unit();
    G4double b = plane.GetB() / plane.GetW().mag();
    Eigen::Vector3d wEigen;
    wEigen << w.x(), w.y(), w.z();
    Eigen::Vector3d EinW = E.ldlt().solve(wEigen);
    G4double wEin = wEigen.transpose() * EinW;
    G4double scale = (b + w.dot(c)) / wEin;
    G4ThreeVector translation = c - scale * G4ThreeVector(EinW(0), EinW(1), EinW(2));
    return translation;
}