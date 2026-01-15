#include "DNAHandler.hh"

#include <vector>
#include <limits>
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAHandler::GetMovedDNA(G4int dnaId, 
                                     const DNAStructure& dnaStructure,
                                     const G4ThreeVector& offset) {
    /* Move the entire DNA structure by a given offset */
    DNA dna = dnaStructure.GetDNA();
    DNAStructure movedDnaStructure;
    
    if(offset==G4ThreeVector(0,0,0)) {
        dna.SetId(dnaId);
        movedDnaStructure.SetDNA(dna);
        return movedDnaStructure;         // return original dna object
    }

    DNA movedDna = DNA(dnaId);

    for (auto& strandPiar : dna.GetStrands()) {
        G4String chainId = strandPiar.first;
        Strand strand = strandPiar.second;
        // make new strand with moved nucleotides
        Strand movedStrand = Strand(chainId);
        for (auto& nucleotidePair : strand.GetNucleotides()) {
            G4int nucleotideId = nucleotidePair.first;
            Nucleotide nucleotide = nucleotidePair.second;
            // make new nucleotide with moved compounds
            Nucleotide movedNucleotide = Nucleotide(nucleotideId, nucleotide.GetBaseType());
            for (auto& compoundPair : nucleotide.GetCompounds()) {
                G4String compoundName = compoundPair.first;
                Compound compound = compoundPair.second;
                // move ellipsoid center
                Ellipsoid ellipsoid = dnaStructure.GetEllipsoid(chainId, nucleotideId, compoundName);
                G4ThreeVector center = ellipsoid.GetCenter();
                ellipsoid.SetCenter(center + offset);
                movedDnaStructure.SetEllipsoid(chainId, nucleotideId, compoundName, ellipsoid);
                // move planes
                std::vector<Plane> planes = dnaStructure.GetPlanes(chainId, nucleotideId, compoundName);
                if (planes.size() != 0) {
                    std::vector<Plane> newPlanes;
                    for (Plane plane : planes) {
                        /* <w,x> + b = 0 → <w,x-a> + b = 0
                           ∴ wx + (b - <w,a>) = 0
                           new b = b - <w,a> */
                        G4double new_b = plane.GetB() - offset.dot(plane.GetW());
                        plane.SetB(new_b);
                        movedDnaStructure.AddPlane(chainId, nucleotideId, compoundName, plane);
                    }
                }
                Compound movedCompound = Compound(compoundName);
                // move atoms
                std::vector<Atom> atoms = compound.GetAtoms();
                for (Atom atom : atoms) {
                    G4ThreeVector atomPos = atom.GetPosition();
                    atom.SetPosition(atomPos + offset);
                    movedCompound.AddAtom(atom);
                }
                // set moved compound to nucleotide
                movedNucleotide.AddCompound(movedCompound);
            }
            // set moved nucleotide to strand
            movedStrand.AddNucleotide(movedNucleotide);
        }
        // set moved strand to dna
        movedDna.AddStrand(movedStrand);
    }
    movedDnaStructure.SetDNA(movedDna);
    return movedDnaStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::pair<G4ThreeVector, G4ThreeVector> DNAHandler::GetRange(const DNA& dna) {
    const G4double kInfinity = std::numeric_limits<G4double>::max();
    G4ThreeVector minPos(kInfinity, kInfinity, kInfinity);
    G4ThreeVector maxPos(-kInfinity, -kInfinity, -kInfinity);

    for (auto& [chainId, strand] : dna.GetStrands()) {
        for (auto& [nucleotideId, nucleotide] : strand.GetNucleotides()) {
            for (auto& [compoundName, compound] : nucleotide.GetCompounds()) {
                for (Atom atom : compound.GetAtoms()) {
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

G4ThreeVector DNAHandler::GetCenter(const DNA& dna) {
    const auto [minPos, maxPos] = GetRange(dna);
    return (minPos + maxPos) * 0.5;
}

void DNAHandler::Print(const DNAStructure& dnaStructure) {
    // dna range calculation
    const DNA dna = dnaStructure.GetDNA();
    const auto [minPos, maxPos] = GetRange(dna);

    G4cout << "==================== DNA " << dna.GetId() << " ====================" << G4endl;
    G4cout << "The range of the dna [Å] : "
           << "[" << minPos.x()/angstrom << ", "
                  << minPos.y()/angstrom << ", "
                  << minPos.z()/angstrom << "], "
           << "[" << maxPos.x()/angstrom << ", "
                  << maxPos.y()/angstrom << ", "
                  << maxPos.z()/angstrom << "]" << G4endl;
    
    for (auto& strandPair : dna.GetStrands()) {
        G4String chainId = strandPair.first;
        Strand strand = strandPair.second;
        G4cout << "================== Strand " << chainId << " ==================" << G4endl;
        for(auto& nucleotidePair : strand.GetNucleotides()) {
            G4int nucleotideId = nucleotidePair.first;
            Nucleotide nucleotide = nucleotidePair.second;
            G4cout << "================== nucleotide " << nucleotideId << " ==================" << G4endl;
            for(auto& compoundPair : nucleotide.GetCompounds()) {
                G4String compoundName = compoundPair.first;
                Compound compound = compoundPair.second;
                std::vector<Atom> atoms = compound.GetAtoms();
                G4String compoundInfo = "Compound: " + compoundName + ", Atoms: ";
                for(size_t i = 0; i < atoms.size(); ++i) {
                    compoundInfo += atoms[i].GetName() + " ";
                }
            G4cout << compoundInfo << G4endl;
            }
        }
    }
    G4cout << "===================================================================" << G4endl;
}