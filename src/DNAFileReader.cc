#include "DNAFileReader.hh"

#include "G4SystemOfUnits.hh"

#include "DNAHandler.hh"

#include <fstream>
#include <iostream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAFileReader::ReadDNAFile(G4int modelId)
{
  /* Read a json file and create a DNAStructure object */
  std::ifstream file(this->fFileName);
  if (!file.is_open())
  {
    G4cerr << "couldn't open json: " << fFileName << G4endl;
    exit(1);
  }
  /* data file format
  {
    "chains": [
      {
        "chain_id": (string),
        "residues": [
          {
            "id": (int),
            "name": (string),
            "compounds": [
              {
                "name": (string),
                "ellipsoid": {
                  "E": (3x3 matrix),
                  "center": [x, y, z]
                },
                "planes": [
                  {
                    "w": [wx, wy, wz],
                    "b": (double)
                  }, ...
                ]
                "atoms": [
                  {
                    "id": (string),
                    "name": (string),
                    "element": (string),
                    "position": [x, y, z],
                    "radius": (double)
                  }, ...
                ]
              }, ...
            ]
          }, ...
        ]
      }, ...
    ]
  }
  */
  json j;
  file >> j;  // parse json file
  file.close();

  // DNA molecule construction
  Model model = GenerateModel(modelId, j);
  // ellipsoid and planes construction
  DNAStructure structure = GenerateStructure(j);
  structure.SetModel(model);

  G4ThreeVector dnaCenter = DNAHandler::GetCenter(model);
  // Generate model object that is moved to the center of gravity
  return DNAHandler::GetMovedStructure(modelId, structure, -dnaCenter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Model DNAFileReader::GenerateModel(G4int modelId, const json& dnaInfo)
{
  Model model = Model(modelId);
  auto chains = dnaInfo["chains"];
  for (const auto& chainInfo : chains)
  {
    // Extract chain ID and create a chain
    G4String chainId = chainInfo["id"].get<G4String>();
    Chain chain = GenerateChain(chainId, chainInfo);
    model.AddChain(chain);
  }
  return model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Chain DNAFileReader::GenerateChain(G4String chainId, const json& chainInfo)
{
  Chain chain = Chain(chainId);
  auto residues = chainInfo["residues"];
  for (const auto& residueInfo : residues)
  {
    // Extract residue ID and create a residue
    G4int residueId = residueInfo["id"];
    Residue residue = GenerateResidue(residueId, residueInfo);
    chain.AddResidue(residue);
  }
  return chain;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue DNAFileReader::GenerateResidue(G4int residueId, const json& residueInfo)
{
  G4String residue_name = residueInfo["name"].get<G4String>();
  Residue residue = Residue(residueId, residue_name);
  auto compounds = residueInfo["compounds"];
  for (const auto& compoundInfo : compounds)
  {
    // Extract compound name and create a compound
    G4String compoundName = compoundInfo["name"].get<G4String>();
    Compound compound = GenerateCompound(compoundName, compoundInfo);
    residue.AddCompound(compound);
  }
  return residue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Compound DNAFileReader::GenerateCompound(G4String compoundName, const json& compoundInfo)
{
  Compound compound = Compound(compoundName);
  json atomsInfo = compoundInfo["atoms"];
  for (const auto& atomInfo : atomsInfo)
  {
    // Create an atom
    Atom atom = GenerateAtom(atomInfo);
    compound.AddAtom(atom);
  }
  return compound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Atom DNAFileReader::GenerateAtom(const json& atomInfo)
{
  G4int id = atomInfo["id"];
  G4String name = atomInfo["name"];
  G4String element = atomInfo["element"];
  std::vector<double> atom_position = atomInfo["position"].get<std::vector<double>>();
  G4double atom_radius = atomInfo["radius"];
  Atom atom = Atom(name, id, element,
                   G4ThreeVector(atom_position[0], atom_position[1], atom_position[2]) * angstrom,
                   atom_radius * angstrom);
  return atom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAFileReader::GenerateStructure(const json& dnaInfo)
{
  /* Parse DNA geometry information such as the oval shape of compounds
     and separating planes from JSON */
  DNAStructure structure = DNAStructure();
  auto chains = dnaInfo["chains"];
  for (const auto& chainInfo : chains)
  {
    G4String chainId = chainInfo["id"].get<G4String>();
    auto residues = chainInfo["residues"];
    for (const auto& residueInfo : residues)
    {
      G4int residueId = residueInfo["id"];
      auto compounds = residueInfo["compounds"];
      for (const auto& compoundInfo : compounds)
      {
        G4String compoundName = compoundInfo["name"].get<G4String>();
        // ellipsoid
        json ellipsoidInfo = compoundInfo["ellipsoid"];
        Ellipsoid ellipsoid = GenerateEllipsoid(ellipsoidInfo);
        structure.SetEllipsoid(chainId, residueId, compoundName, ellipsoid);

        // planes
        json planesInfo = compoundInfo["planes"];
        for (auto& planeInfo : planesInfo)
        {
          Plane plane = GeneratePlane(planeInfo);
          structure.AddPlane(chainId, residueId, compoundName, plane);
        }
      }
    }
  }
  return structure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid DNAFileReader::GenerateEllipsoid(const json& ellipsoidInfo)
{
  std::vector<std::vector<G4double>> E_matrix =
    ellipsoidInfo["E"].get<std::vector<std::vector<G4double>>>();
  std::vector<G4double> center = ellipsoidInfo["center"].get<std::vector<G4double>>();

  // convert to matrix using Eigen
  Eigen::Matrix3d mat;
  for (G4int i = 0; i < 3; ++i)
    for (G4int j = 0; j < 3; ++j)
      mat(i, j) = E_matrix[i][j];


  // center of the ellipsoid
  G4ThreeVector centerVector = G4ThreeVector(center[0], center[1], center[2]) * angstrom;

  return Ellipsoid(mat, centerVector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Plane DNAFileReader::GeneratePlane(const json& planeInfo)
{
  std::vector<G4double> w = planeInfo["w"].get<std::vector<G4double>>();
  G4double b = planeInfo["b"];
  return Plane(G4ThreeVector(w[0], w[1], w[2]), b * angstrom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
