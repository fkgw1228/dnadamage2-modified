#include "DNAFileReader.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include "G4SystemOfUnits.hh"
#include "DNAHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAFileReader::ReadDNAFile(G4int dnaId) {
    /* Read a json file and create a DNAStructure object */
    std::ifstream file(this->fFileName);
    if(!file.is_open()) {
        G4cerr << "couldn't open json: " << fFileName << G4endl;
        exit(1);
    }
    /* data file format 
    {
        "strands": [
            {
                "chain_id": (string),
                "nucleotides": [
                    {
                        "id": (int),
                        "base_type": (string),
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
    file >> j;      // parse json file
    file.close();

    // DNA molecule construction
    DNA dna = GenerateDNA(dnaId, j);
    // ellipsoid and planes construction
    DNAStructure dnaStructure = GenerateGeometry(j);
    dnaStructure.SetDNA(dna);

    G4ThreeVector dnaCenter = DNAHandler::GetCenter(dna);
    // Generate dna object that is moved to the center of gravity
    return DNAHandler::GetMovedDNA(dnaId, dnaStructure, -dnaCenter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNA DNAFileReader::GenerateDNA(G4int dnaId, const json& dnaInfo) {
    DNA dna = DNA(dnaId);
    auto strands = dnaInfo["strands"];
    for (const auto& strandInfo : strands) {
        // Extract chain ID and create a strand
        G4String chainId = strandInfo["chain_id"].get<G4String>();
        Strand strand = GenerateStrand(chainId, strandInfo);
        dna.AddStrand(strand);
    }
    return dna;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Strand DNAFileReader::GenerateStrand(G4String chainId, const json& strandInfo) {
    Strand strand = Strand(chainId);
    auto nucleotides = strandInfo["nucleotides"];
    for (const auto& nucleotideInfo : nucleotides) {
        // Extract nucleotide ID and create a nucleotide
        G4int nucleotideId = nucleotideInfo["id"];
        Nucleotide nucleotide = GenerateNucleotide(nucleotideId, nucleotideInfo);
        strand.AddNucleotide(nucleotide);
    }
    return strand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Nucleotide DNAFileReader::GenerateNucleotide(G4int nucleotideId, const json& nucleotideInfo) {
    G4String base_type = nucleotideInfo["base_type"].get<G4String>();
    Nucleotide nucleotide = Nucleotide(nucleotideId, base_type);
    auto compounds = nucleotideInfo["compounds"];
    for (const auto& compoundInfo : compounds) {
        // Extract compound name and create a compound
        G4String compoundName = compoundInfo["name"].get<G4String>();
        Compound compound = GenerateCompound(compoundName, compoundInfo);
        nucleotide.AddCompound(compound);
    }
    return nucleotide;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Compound DNAFileReader::GenerateCompound(G4String compoundName, const json& compoundInfo) {
    Compound compound = Compound(compoundName);
    json atomsInfo = compoundInfo["atoms"];
    for (const auto& atomInfo : atomsInfo) {
        // Create an atom
        Atom atom = GenerateAtom(atomInfo);
        compound.AddAtom(atom);
    }
    return compound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Atom DNAFileReader::GenerateAtom(const json& atomInfo) {
    G4int id = atomInfo["id"];
    G4String name = atomInfo["name"];
    G4String element = atomInfo["element"];
    std::vector<double> atom_position = atomInfo["position"].get<std::vector<double>>();
    G4double atom_radius = atomInfo["radius"];
    Atom atom = Atom(name, id, element,
                     G4ThreeVector(atom_position[0],
                                   atom_position[1],
                                   atom_position[2]) * angstrom,
                     atom_radius * angstrom);
    return atom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAStructure DNAFileReader::GenerateGeometry(const json& dnaInfo) {
    /* Parse DNA geometry information such as the oval shape of compounds
       and separating planes from JSON */
    DNAStructure dnaStructure = DNAStructure();
    auto strands = dnaInfo["strands"];
    for (const auto& strandInfo : strands) {
        G4String chainId = strandInfo["chain_id"].get<G4String>();
        auto nucleotides = strandInfo["nucleotides"];
        for (const auto& nucleotideInfo : nucleotides) {
            G4int nucleotideId = nucleotideInfo["id"];
            auto compounds = nucleotideInfo["compounds"];
            for (const auto& compoundInfo : compounds) {
                G4String compoundName = compoundInfo["name"].get<G4String>();
                // ellipsoid
                json ellipsoidInfo = compoundInfo["ellipsoid"];
                Ellipsoid ellipsoid = GenerateEllipsoid(ellipsoidInfo);
                dnaStructure.SetEllipsoid(chainId, nucleotideId, compoundName, ellipsoid);

                // planes
                json planesInfo = compoundInfo["planes"];
                for (auto& planeInfo : planesInfo) {
                    Plane plane = GeneratePlane(planeInfo);
                    dnaStructure.AddPlane(chainId, nucleotideId, compoundName, plane);
                }
            }
        }
    }
    return dnaStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Ellipsoid DNAFileReader::GenerateEllipsoid(const json& ellipsoidInfo) {
    
    std::vector<std::vector<G4double>> E_matrix = ellipsoidInfo["E"].get<std::vector<std::vector<G4double>>>();
    std::vector<G4double> center = ellipsoidInfo["center"].get<std::vector<G4double>>();

    // convert to matrix using Eigen
    Eigen::Matrix3d mat;
    mat << E_matrix[0][0], E_matrix[0][1], E_matrix[0][2],
           E_matrix[0][1], E_matrix[1][1], E_matrix[1][2],
           E_matrix[0][2], E_matrix[1][2], E_matrix[2][2];
    
    // center of the ellipsoid
    G4ThreeVector centerVector = G4ThreeVector(center[0], center[1], center[2]) * angstrom;

    return Ellipsoid(mat, centerVector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Plane DNAFileReader::GeneratePlane(const json& planeInfo) {
    std::vector<G4double> w = planeInfo["w"].get<std::vector<G4double>>();
    G4double b = planeInfo["b"];
    return Plane(G4ThreeVector(w[0], w[1], w[2]), b * angstrom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......