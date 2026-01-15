#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "nlohmann/json.hpp"
#include "Eigen/Dense"

#include "TCanvas.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TH1F.h"

using json = nlohmann::json;

void PlotDamageLocation(std::string jsonFileName,
                        std::string ssbDsbFileName,
                        std::string outputImageFileName);
void LoadJsonFile(std::string jsonFileName);
void LoadDamageFile(std::string strandBreakFileName);
void DrawDNAMolecule(TPad* canvas, int axis);
void DrawDamageLocations(TPad* canvas, int axis);
Eigen::Matrix2d Get2DEllipsoidMatrix(Eigen::Matrix3d E, int axis);
double GetRotationTheta(Eigen::Matrix2d E2D);
Eigen::Vector2d Get2DPosition(Eigen::Vector3d vec, int axis);

struct DNALocation {
    std::string chainId;
    int nucleotideId;
    std::string compoundName;

    bool operator<(const DNALocation& other) const {
        if (chainId != other.chainId) return chainId < other.chainId;
        if (nucleotideId != other.nucleotideId) return nucleotideId < other.nucleotideId;
        return compoundName < other.compoundName;
    }
};

struct AtomInfo {
    std::string atomName;
    Eigen::Vector3d position;
};

struct Ellipsoid {
    Eigen::Matrix3d E;
    Eigen::Vector3d center;
};

struct DamageInfo {
    std::string atomOrMoleculeName;
    int breakId;
};

std::map<DNALocation, std::vector<AtomInfo>> atomMap;
std::map<DNALocation, Ellipsoid> ellipsoidMap;
std::map<DNALocation, std::vector<DamageInfo>> damageAtomMap;
double minX = 1e9, maxX = -1e9;
double minY = 1e9, maxY = -1e9;
double minZ = 1e9, maxZ = -1e9;

int main(int argc, char** argv) {
    if (argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <json_file_name> <ssb_dsb_file_name> [output_image_file_name]"
                  << std::endl;
        return 1;
    }
    std::string jsonFileName = argv[1];
    std::string ssbDsbFileName = argv[2];
    std::string outputImageFileName = argc > 3 ? argv[3] : "";

    try {
        PlotDamageLocation(jsonFileName, ssbDsbFileName, outputImageFileName);
    } catch (const std::exception& e) {
        std::cerr << "Error during plotting: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

void PlotDamageLocation(std::string jsonFileName,
                        std::string strandBreakFileName,
                        std::string outputImageFileName) {
    // canvas
    TCanvas* canvas = new TCanvas("canvas","Damage Locations",800,800);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetTopMargin(0.15);
    canvas->SetBottomMargin(0.15);

    int axis = 1; // 0=x, yz plane; 1=y, xz plane; 2=z, xy plane

    LoadJsonFile(jsonFileName);
    LoadDamageFile(strandBreakFileName);

    // Set axis labels based on projection plane
    const char* xLabel = "";
    const char* yLabel = "";
    if (axis == 0) {            // 0=x, yz plane
        xLabel = "Y [#AA]";
        yLabel = "Z [#AA]";
    } else if (axis == 1) {     // 1=y, xz plane
        xLabel = "X [#AA]";
        yLabel = "Z [#AA]";
    } else if (axis == 2) {     // 2=z, xy plane
        xLabel = "X [#AA]";
        yLabel = "Y [#AA]";
    }
    
    // Draw frame with axis labels on canvas
    TH1F* frame = canvas->DrawFrame(-15, -15, 15, 15);
    frame->GetXaxis()->SetTitle(xLabel);
    frame->GetYaxis()->SetTitle(yLabel);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    
    // Create pad on top for drawing
    TPad* plotPad = new TPad("plotPad", "", 0.15, 0.15, 0.85, 0.85);
    plotPad->SetFillStyle(0);  // Transparent
    plotPad->SetBit(TPad::kClipFrame);
    plotPad->Draw();
    plotPad->cd();
    plotPad->Range(-15, -15, 15, 15);
    plotPad->SetGrid();
    
    DrawDNAMolecule(plotPad, axis);
    DrawDamageLocations(plotPad, axis);

    canvas->cd();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(outputImageFileName.c_str());
    delete canvas;  // Clean up memory
    return;
}

void LoadJsonFile(std::string jsonFileName) {
    std::ifstream jsonFile(jsonFileName);
    if (!jsonFile.is_open()) throw std::runtime_error("Could not open JSON file " + jsonFileName);
    json j;
    jsonFile >> j;
    jsonFile.close();

    auto& strands = j["strands"];
    for (auto& strand : strands) {
        std::string chainId = strand["chain_id"];
        auto& nucleotides = strand["nucleotides"];
        for (auto& nucleotide : nucleotides) {
            int nucleotideId = nucleotide["id"];
            auto& compounds = nucleotide["compounds"];
            for (auto& compound : compounds) {
                std::string compoundName = compound["name"];
                DNALocation location{chainId, nucleotideId, compoundName};

                auto& ell = compound["ellipsoid"];
                Eigen::Matrix3d E;
                for (int i=0; i<3; i++) for (int j=0; j<3; j++) E(i,j) = ell["E"][i][j];
                Eigen::Vector3d center;
                for (int i=0; i<3; i++) center(i) = ell["center"][i];
                ellipsoidMap[location] = Ellipsoid{E, center};
                
                for (auto& atom : compound["atoms"]) {
                    auto pos = atom["position"];
                    AtomInfo ai{atom["name"], Eigen::Vector3d(pos[0], pos[1], pos[2])};
                    atomMap[location].push_back(ai);
                    minX = std::min(minX, ai.position(0)); maxX = std::max(maxX, ai.position(0));
                    minY = std::min(minY, ai.position(1)); maxY = std::max(maxY, ai.position(1));
                    minZ = std::min(minZ, ai.position(2)); maxZ = std::max(maxZ, ai.position(2));
                }
            }
        }
    }
    // Move center of DNA to origin
    Eigen::Vector3d dnaCenter = Eigen::Vector3d((minX+maxX)*0.5, (minY+maxY)*0.5, (minZ+maxZ)*0.5);
    for (auto& [location, ellipsoid] : ellipsoidMap)
        ellipsoid.center -= dnaCenter;
    for (auto& [location, atomList] : atomMap)
        for (auto& atom : atomList)
            atom.position -= dnaCenter;
    return;
}

void LoadDamageFile(std::string strandBreakFileName) {
    std::ifstream damageFile(strandBreakFileName);
    if (!damageFile.is_open()) throw std::runtime_error("Could not open damage file " + strandBreakFileName);
    std::string line;
    getline(damageFile, line); getline(damageFile, line); // skip header
    while (std::getline(damageFile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string chainId, compoundName, atomOrMoleculeName;
        int eventId, dnaId, nucleotideId, breakId;
        double dose;
        if (!(iss >> eventId >> dnaId >> chainId >> nucleotideId
                  >> compoundName >> atomOrMoleculeName >> breakId >> dose)) {
            continue; // skip malformed lines
        }
        if (atomOrMoleculeName == "H1") atomOrMoleculeName = "H1'";
        if (atomOrMoleculeName == "H2a") atomOrMoleculeName = "H2'";
        if (atomOrMoleculeName == "H2b") atomOrMoleculeName = "H2''";
        if (atomOrMoleculeName == "H3") atomOrMoleculeName = "H3'";
        if (atomOrMoleculeName == "H4") atomOrMoleculeName = "H4'";
        if (atomOrMoleculeName == "H5a") atomOrMoleculeName = "H5'";
        if (atomOrMoleculeName == "H5b") atomOrMoleculeName = "H5''";
        if (dnaId != 1) continue; // only process dnaId 1
        DNALocation location{chainId, nucleotideId, compoundName};
        DamageInfo damageInfo = {atomOrMoleculeName, breakId};
        damageAtomMap[location].push_back(damageInfo);
    }
    return;
}

void DrawDNAMolecule(TPad* canvas, int axis) {
    for (const auto& [location, ellipsoid] : ellipsoidMap) {
        Eigen::Matrix3d E = ellipsoid.E;
        Eigen::Vector3d center = ellipsoid.center;

        Eigen::Matrix2d E2D = Get2DEllipsoidMatrix(E, axis);
        Eigen::Vector2d center2D = Get2DPosition(center, axis);
        double a = 1.0 / sqrt(E2D(0,0));
        double b = 1.0 / sqrt(E2D(1,1));
        double theta = GetRotationTheta(E2D);
        
        TEllipse* ellipse = new TEllipse(center2D(0), center2D(1), a, b, 0, 360, theta);
        ellipse->SetFillStyle(1);
        ellipse->SetLineWidth(1);
        int color = kBlack;
        if (location.compoundName == "phosphate")  color = kYellow+1;
        else if (location.compoundName == "deoxyribose") color = kRed-7;
        else if (location.compoundName == "base")  color = kGreen-6;
        ellipse->SetFillColorAlpha(color, 0.2);
        ellipse->SetLineColor(color);
        ellipse->Draw("same");
    }
    return;
}

void DrawDamageLocations(TPad* canvas, int axis) {
    for (const auto& [location, damagedAtoms] : damageAtomMap) {
        if (atomMap.find(location) == atomMap.end()) continue;
        const auto& atoms = atomMap[location];
        for (const auto& atomInfo : atoms) {
            for (const auto& damagedAtomInfo : damagedAtoms) {
                if (atomInfo.atomName != damagedAtomInfo.atomOrMoleculeName) continue;
                Eigen::Vector3d position3D = atomInfo.position;
                Eigen::Vector2d position2D = Get2DPosition(position3D, axis);

                TGraph* point = new TGraph(1);
                point->SetPoint(0, position2D(0), position2D(1));
                int color = kBlack;
                if (damagedAtomInfo.breakId == 1) color = kRed;
                else if (damagedAtomInfo.breakId == 2) color = kBlue;
                point->SetMarkerStyle(20);
                point->SetMarkerSize(1.5);
                point->SetMarkerColor(color);
                point->Draw("P same");
            }
        }
    }
    return;
}

Eigen::Matrix2d Get2DEllipsoidMatrix(Eigen::Matrix3d E, int axis) {
    // axis: 0=x, 1=y, 2=z
    Eigen::Matrix2d E2D;
    if (axis == 0) {               // yz plane
        E2D(0,0) = E(1,1); E2D(0,1) = E(1,2);
        E2D(1,0) = E(2,1); E2D(1,1) = E(2,2);
    } else if (axis == 1) {        // xz plane
        E2D(0,0) = E(0,0); E2D(0,1) = E(0,2);
        E2D(1,0) = E(2,0); E2D(1,1) = E(2,2);
    } else if (axis == 2) {        // xy plane
        E2D(0,0) = E(0,0); E2D(0,1) = E(0,1);
        E2D(1,0) = E(1,0); E2D(1,1) = E(1,1);
    }
    return E2D;
}

double GetRotationTheta(Eigen::Matrix2d E2D) {
    /* Rotation angle theta is defilned as:
        theta = 0.5 * arctan(2b / (a - c))
       where the ellipsoid matrix is given by
        | a  b |
        | b  c |
    */
    double a = E2D(0,0);
    double b = E2D(0,1);
    double c = E2D(1,1);
    double theta = 0.5 * atan2(2.0 * b, a - c);
    return theta * 180.0 / M_PI + 90.0; // convert to degrees
}

Eigen::Vector2d Get2DPosition(Eigen::Vector3d vec, int axis) {
    // axis: 0=x, 1=y, 2=z
    Eigen::Vector2d vec2D;
    if (axis == 0) {               // yz plane
        vec2D(0) = vec(1);
        vec2D(1) = vec(2);
    } else if (axis == 1) {        // xz plane
        vec2D(0) = vec(0);
        vec2D(1) = vec(2);
    } else if (axis == 2) {        // xy plane
        vec2D(0) = vec(0);
        vec2D(1) = vec(1);
    }
    return vec2D;
}