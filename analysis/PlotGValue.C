#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TH1F.h>

void PlotGValue(TString filename) {
    if (filename == "") {
        printf("Please provide a valid .root filename.");
        printf("Usage Example: .x PlotGValue(\"data/species/SpeciesInfo.root\")");
        return;
    }

    // No statistics box
    gStyle->SetOptStat(0);

    // Open the ROOT file
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        printf("Error: Could not open file %s\n", filename.Data());
        return;
    }

    // Get TTree
    TTree *tree = (TTree*)file->Get("species");
    if (!tree) {
        printf("Error: tree not found in file %s\n", filename.Data());
        return;
    }

    // Bind variables
    Double_t time;
    Double_t G;
    Double_t RMS;
    Double_t Nb;
    char     Molecule[256];

    tree->SetBranchAddress("Time", &time);
    tree->SetBranchAddress("G", &G);
    tree->SetBranchAddress("RMS", &RMS);
    tree->SetBranchAddress("Nb", &Nb);
    tree->SetBranchAddress("Molecule", Molecule);

    std::vector<TString> speciesNames;
    for (Long64_t i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        TString speciesNameStr(Molecule);
        if (std::find(speciesNames.begin(), speciesNames.end(), speciesNameStr) == speciesNames.end()) {
            speciesNames.push_back(speciesNameStr);
        }
    }
    printf("Species in the tree:\n");
    for (const auto& sname : speciesNames) {
        printf(" - %s\n", sname.Data());
    }

    std::vector<TString> plotSpecies1 = {
        "e_aq^-1",
        "°OH^0",
        "H^0",
        "H2O2^0",
        "HO_2°^0",
        "O_2^-1"
    };
    std::vector<Int_t> plotMarkers1 = {
        24, 20, 26, 22, 25, 21
    };
    std::vector<TString> legendNames1 = {
        "e_aq",
        "OH#bullet",
        "H#bullet",
        "H_{2}O_{2}",
        "HO_{2}^{#bullet}",
        "O_{2}^{#bullet-}"
    };
    std::vector<TString> plotSpecies2 = {
        "DNA_DamagedDeoxyribose_H1_OH^0",
        "DNA_DamagedDeoxyribose_H2a_OH^0",
        "DNA_DamagedDeoxyribose_H2b_OH^0",
        "DNA_DamagedDeoxyribose_H3_OH^0",
        "DNA_DamagedDeoxyribose_H4_OH^0",
        "DNA_DamagedDeoxyribose_H5a_OH^0",
        "DNA_DamagedDeoxyribose_H5b_OH^0"
    };
    std::vector<Int_t> plotMarkers2 = {
        20, 24, 26, 22, 25, 21, 30
    };
    std::vector<TString> legendNames2 = {
        "H_{1}",
        "H_{2a}",
        "H_{2b}",
        "H_{3}",
        "H_{4}",
        "H_{5a}",
        "H_{5b}"
    };

    // Create canvas
    TMultiGraph *mg1 = new TMultiGraph();
    TCanvas *c1 = new TCanvas("c1", "GValue vs Time", 800, 800);
    auto leg1 = new TLegend(0.7, 0.6, 0.9, 0.9, "#font[22]{chemical species}");
    leg1->SetBorderSize(0); // No border
    leg1->SetFillStyle(0);  // Transparent background

    gPad->SetLogx(true);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.07); 
    for (Int_t i=0; i<plotSpecies1.size(); i++) {
        TString species = plotSpecies1[i];
        TGraphErrors *gr = new TGraphErrors();
        for (Long64_t j=0; j<tree->GetEntries(); j++) {
            tree->GetEntry(j);
            if (!std::isfinite(G)) {
                continue;
            }
            TString speciesNameStr(Molecule);
            if (speciesNameStr == species) {
                size_t index = gr->GetN();
                gr->SetPoint(index, time, G);
                gr->SetPointError(index, 1e-10, RMS);
            }
        }
        gr->SetMarkerStyle(plotMarkers1[i]);
        gr->SetLineStyle(1);
        gr->SetMarkerSize(0.5);
        gr->SetLineWidth(1);
        mg1->Add(gr, "LPE");
        leg1->AddEntry(gr, legendNames1[i].Data(), "lpe");
    }
    mg1->SetTitle(";Time [ps]; G Value [molecules/100 eV]");
    mg1->Draw("APLE");
    leg1->Draw();
    gPad->Update();
    mg1->GetXaxis()->SetRangeUser(1, 1e7);
    mg1->GetXaxis()->SetTitleOffset(1.5);
    mg1->GetYaxis()->SetTitleOffset(1.5);
    c1->Modified();
    c1->Update();

    std::string res;
    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        c1->SaveAs(outFilename.Data());
    }

    TMultiGraph *mg2 = new TMultiGraph();
    TCanvas *c2 = new TCanvas("c2", "GValue vs Time", 800, 800);
    auto leg2 = new TLegend(0.2, 0.6, 0.4, 0.9,
                            "#font[22]{damaged species}");
    leg2->SetBorderSize(0); // No border
    leg2->SetFillStyle(0);  // Transparent background

    gPad->SetLogx(true);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.07); 
    for (Int_t i=0; i<plotSpecies2.size(); i++) {
        TString species = plotSpecies2[i];
        TGraphErrors *gr = new TGraphErrors();
        gr->SetTitle(";Time [ps]; G Value [molecules/100 eV]");
        for (Long64_t j=0; j<tree->GetEntries(); j++) {
            tree->GetEntry(j);
            TString speciesNameStr(Molecule);
            if (speciesNameStr == species) {
                size_t index = gr->GetN();
                gr->SetPoint(index, time, G);
                gr->SetPointError(index, 1e-10, RMS);
            }
        }
        gr->SetMarkerStyle(plotMarkers2[i]);
        gr->SetLineStyle(1);
        gr->SetMarkerSize(0.5);
        gr->SetLineWidth(1);
        mg2->Add(gr, "LPE");
        leg2->AddEntry(gr, legendNames2[i].Data(), "lpe");
    }
    mg2->SetTitle(";Time [ps]; G Value [molecules/100 eV]");
    mg2->Draw("APLE");
    leg2->Draw();
    gPad->Update();
    mg2->GetXaxis()->SetRangeUser(1, 1e7);
    mg2->GetXaxis()->SetTitleOffset(1.5);
    mg2->GetYaxis()->SetTitleOffset(2);
    c2->Modified();
    c2->Update();

    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        c1->SaveAs(outFilename.Data());
    }
}