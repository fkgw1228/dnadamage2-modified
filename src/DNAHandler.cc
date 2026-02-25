#include "DNAHandler.hh"

#include "G4SystemOfUnits.hh"

#include <limits>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAHandler::GetMovedStructure(G4int structureId,
                                        DNAStructure& structure,
                                        G4ThreeVector offset)
{
  /* Move the entire DNA structure by a given offset */
  Model model = structure.GetModel();
  DNAStructure movedDnaStructure;

  Model movedModel = Model(structureId);

  for (auto& chainPair : model.GetChains())
  {
    G4String chainId = chainPair.first;
    Chain chain = chainPair.second;
    // make new chain with moved residues
    Chain movedChain = Chain(chainId);
    for (auto& residuePair : chain.GetResidues())
    {
      G4int residueId = residuePair.first;
      Residue residue = residuePair.second;
      // make new residue with moved compounds
      Residue movedResidue = Residue(residueId, residue.GetName());
      for (auto& compoundPair : residue.GetCompounds())
      {
        G4String compoundName = compoundPair.first;
        Compound compound = compoundPair.second;
        // move ellipsoid center
        Ellipsoid ellipsoid = structure.GetEllipsoid(chainId, residueId, compoundName);
        G4ThreeVector center = ellipsoid.GetCenter();
        ellipsoid.SetCenter(center + offset);
        movedDnaStructure.SetEllipsoid(chainId, residueId, compoundName, ellipsoid);
        // move planes
        std::vector<Plane> planes = structure.GetPlanes(chainId, residueId, compoundName);
        if (planes.size() != 0)
        {
          std::vector<Plane> newPlanes;
          for (Plane plane : planes)
          {
            /* <w,x> + b = 0 → <w,x-a> + b = 0
               ∴ wx + (b - <w,a>) = 0
               new b = b - <w,a> */
            G4double new_b = plane.GetB() - offset.dot(plane.GetW());
            plane.SetB(new_b);
            movedDnaStructure.AddPlane(chainId, residueId, compoundName, plane);
          }
        }
        Compound movedCompound = Compound(compoundName);
        // move atoms
        std::vector<Atom> atoms = compound.GetAtoms();
        for (Atom atom : atoms)
        {
          G4ThreeVector atomPos = atom.GetPosition();
          atom.SetPosition(atomPos + offset);
          movedCompound.AddAtom(atom);
        }
        // set moved compound to residue
        movedResidue.AddCompound(movedCompound);
      }
      // set moved residue to chain
      movedChain.AddResidue(movedResidue);
    }
    // set moved chain to model
    movedModel.AddChain(movedChain);
  }
  movedDnaStructure.SetModel(movedModel);
  return movedDnaStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::pair<G4ThreeVector, G4ThreeVector> DNAHandler::GetRange(Model& model)
{
  const G4double kInfinity = std::numeric_limits<G4double>::max();
  G4ThreeVector minPos(kInfinity, kInfinity, kInfinity);
  G4ThreeVector maxPos(-kInfinity, -kInfinity, -kInfinity);

  for (auto& [chainId, chain] : model.GetChains())
  {
    for (auto& [residueId, residue] : chain.GetResidues())
    {
      for (auto& [compoundName, compound] : residue.GetCompounds())
      {
        for (Atom atom : compound.GetAtoms())
        {
          G4ThreeVector pos = atom.GetPosition();
          minPos.setX(std::min(minPos.x(), pos.x()));
          minPos.setY(std::min(minPos.y(), pos.y()));
          minPos.setZ(std::min(minPos.z(), pos.z()));
          maxPos.setX(std::max(maxPos.x(), pos.x()));
          maxPos.setY(std::max(maxPos.y(), pos.y()));
          maxPos.setZ(std::max(maxPos.z(), pos.z()));
        }
      }
    }
  }

  return {minPos, maxPos};
}

G4ThreeVector DNAHandler::GetCenter(Model& dna)
{
  const auto [minPos, maxPos] = GetRange(dna);
  return (minPos + maxPos) * 0.5;
}

void DNAHandler::Print(DNAStructure& structure)
{
  // dna range calculation
  Model model = structure.GetModel();
  auto [minPos, maxPos] = GetRange(model);

  G4cout << "==================== Model " << model.GetId() << " ====================" << G4endl;
  G4cout << "The range of the dna [Å] : "
         << "[" << minPos.x() / angstrom << ", " << minPos.y() / angstrom << ", "
         << minPos.z() / angstrom << "], "
         << "[" << maxPos.x() / angstrom << ", " << maxPos.y() / angstrom << ", "
         << maxPos.z() / angstrom << "]" << G4endl;

  for (auto& chainPair : model.GetChains())
  {
    G4String chainId = chainPair.first;
    Chain chain = chainPair.second;
    G4cout << "================== Chain " << chainId << " ==================" << G4endl;
    for (auto& residuePair : chain.GetResidues())
    {
      G4int residueId = residuePair.first;
      Residue residue = residuePair.second;
      G4cout << "================== residue " << residueId << " ==================" << G4endl;
      for (auto& compoundPair : residue.GetCompounds())
      {
        G4String compoundName = compoundPair.first;
        Compound compound = compoundPair.second;
        std::vector<Atom> atoms = compound.GetAtoms();
        G4String compoundInfo = "Compound: " + compoundName + ", Atoms: ";
        for (size_t i = 0; i < atoms.size(); ++i)
        {
          compoundInfo += atoms[i].GetName() + " ";
        }
        G4cout << compoundInfo << G4endl;
      }
    }
  }
  G4cout << "===================================================================" << G4endl;
}