/*
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/


#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>
#include <TMath.h>
#include <TH1F.h>
#include <vector>
#include <limits>
#include <TString.h>
#include <stdio.h>
#include <TFile.h>
#include <map>
#include <algorithm>

void PlotSpeciesPos(TString filename) {
    if (filename == "") {
        printf("Please provide a valid root filename.\n");
        printf("Usage example: .x PlotSpeciesPos(\"data/species/SpeciesDist.root\")");
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
    TTree *tree = (TTree*)file->Get("SpeciesTimeDistribution");
    if (!tree) {
        printf("Error: tree not found in file %s\n", filename.Data());
        return;
    }

    // Bind variables
    Int_t    eventID;
    Int_t    trackID;
    Double_t time_ps;
    Double_t x_ang, y_ang, z_ang;
    char     name[256];
    tree->SetBranchAddress("EventID", &eventID);
    tree->SetBranchAddress("TrackID", &trackID);
    tree->SetBranchAddress("Time", &time_ps);
    tree->SetBranchAddress("PosX", &x_ang);
    tree->SetBranchAddress("PosY", &y_ang);
    tree->SetBranchAddress("PosZ", &z_ang);
    tree->SetBranchAddress("SpeciesName", name);

    // Print Species Names in the tree
    std::vector<TString> speciesNames;
    for (Long64_t i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        TString speciesNameStr(name);
        if (std::find(speciesNames.begin(), speciesNames.end(), speciesNameStr) == speciesNames.end()) {
            speciesNames.push_back(speciesNameStr);
        }
    }
    printf("Species in the tree:\n");
    for (const auto& sname : speciesNames) {
        printf(" - %s\n", sname.Data());
    }

    // Total entries
    Long64_t nentries = tree->GetEntries();
    if (nentries == 0) {
        printf("Tree is empty in file %s\n", filename.Data());
        return;
    }

    // Calculate and print entries for DNA and other species
    Long64_t nDna = tree->GetEntries("strstr(SpeciesName, \"DNA_Deoxyribose\")  || strstr(SpeciesName, \"DNA_Histone\")");
    Long64_t nSpecies = tree->GetEntries("!strstr(SpeciesName, \"DNA_Deoxyribose\") && !strstr(SpeciesName, \"DNA_Histone\")");

    std::cout << "Total entries in DNA tree: " << nDna << std::endl;
    std::cout << "Total entries in Species tree: " << nSpecies << std::endl;

    std::map<Int_t, TString> prefixMap ={
        {0, "p"},
        {1, "n"},
        {2, "u"},
        {3, "m"}
    };

    std::map<TString, std::vector<Int_t>> speciesEntryIDMap;
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        TString speciesNameStr(name);
        speciesEntryIDMap[speciesNameStr].push_back(i);
    }

    // divide time into bins
    const int nPads = 8;
    const int nPointsPerGraph = 10000 + nDna;

    // Provide 9 edges for 8 pads: [t0,t1), [t1,t2), ... [t7,t8)
    std::vector<Double_t> timeEdges = {
        // 1ps, 10ps, 100ps, 1ns, 5ns, 10ns, 100ns, 1us, 10us
        1.0, 10.0, 100.0, 1000.0, 5000.0, 10000.0, 100000.0, 1000000.0, 10000000.0, 100000000.0
    };


    // Create canvas with multiple pads
    TCanvas *c = new TCanvas("c", "SpeciesDist at each time", 1600, 800);
    c->Divide(4, 2); // 4 columns, 2 rows
    std::vector<TGraph*> graphsDamaged;
    
    for (Int_t padIndex=0; padIndex<nPads; padIndex++) {
        if (padIndex > 7) break; // canvas is 4x2
        c->cd(padIndex+1);

        int numPointsPlotted = 0;

        // Increase margins to prevent axis labels from being cut off
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.1);

        // set the size of the each pad to be square
        gPad->SetFixedAspectRatio();

        Double_t timeLow = timeEdges[padIndex];
        Double_t timeHigh = timeEdges[padIndex+1];

        double_t log10_timeLow = TMath::Log10(timeLow);
        double_t log10_timeHigh = TMath::Log10(timeHigh);
        double_t timeLowForPrint =  (double_t)(TMath::Power(10, log10_timeLow-(Int_t)log10_timeLow)) * (double_t)TMath::Power(10, (int)TMath::Log10(timeLow)%3);
        double_t timeHighForPrint =   (double_t)(TMath::Power(10, log10_timeHigh-(Int_t)log10_timeHigh)) * (double_t)TMath::Power(10, (int)TMath::Log10(timeHigh)%3);
        TString prefix1 = prefixMap[(int)TMath::Log10(timeLow)/3];
        TString prefix2 = prefixMap[(int)TMath::Log10(timeHigh)/3];
        printf("Pad %d : Time range %.2f %ss - %.2f %ss\n", padIndex, timeLowForPrint, prefix1.Data(), timeHighForPrint, prefix2.Data());

        // Create an empty frame first (so axes exist even when a time bin has 0 points)
        TH1F* frame = (TH1F*)gPad->DrawFrame(
            -100, -150, 100, 150,
            Form("Time: %.2f %ss - %.2f %ss",
                 timeLowForPrint, prefix1.Data(), timeHighForPrint, prefix2.Data())
        );
        frame->GetXaxis()->SetTitle("x [#AA]");
        frame->GetYaxis()->SetTitle("z [#AA]");
        frame->GetXaxis()->SetNoExponent();
        frame->GetYaxis()->SetNoExponent();

        std::map<Int_t, std::vector<Int_t>> addedTrackIDsPerEvent;

        std::vector<TGraph*> graphsNormal;

        for (auto [name, entryIDs] : speciesEntryIDMap) {

            TGraph *g1 = new TGraph();
            g1->SetMarkerSize(0.4);
            g1->SetMarkerStyle(20);

            bool isDamaged = name.Contains("DNA_DamagedDeoxyribose");
            bool isDeoxyribose = name.Contains("DNA_Deoxyribose");
            bool isHistone = name.Contains("DNA_Histone");

            if      (name == "°OH^0") {
                g1->SetMarkerColor(kCyan-3);
            }
            else if (isDamaged) {
                g1->SetMarkerColor(kBlack);
            }
            else if (isDeoxyribose) {
                g1->SetMarkerColor(kRed-7);
            }
            else if (isHistone) {
                g1->SetMarkerColor(kOrange-3); // Histone color
            }
            else  {
                g1->SetMarkerColor(kGreen+2);
            }
            
            for (auto entryID : entryIDs) {
                tree->GetEntry(entryID);

                bool shouldPlot = false;

                if ((isDeoxyribose || isHistone) && !isDamaged) {
                    // Static DNA/Histone structure: only plot Event 0 at 1 ps on ALL pads
                    if (eventID == 0 && TMath::Abs(time_ps - 1.0) < 0.001) {
                        shouldPlot = true;
                    }
                } else {
                    // Dynamic species: plot if they exist in the current time window
                    if (time_ps >= timeLow && time_ps < timeHigh) {
                        shouldPlot = true;
                    }
                }

                if (shouldPlot) {
                    //if (std::find(addedTrackIDsPerEvent[eventID].begin(), addedTrackIDsPerEvent[eventID].end(), trackID) != addedTrackIDsPerEvent[eventID].end()) {
                    //    continue; // skip if this trackID is already added in this time bin
                    //}
                    // if (padIndex == 0 && !isDeoxyribose) {
                    //     printf("Adding point: Name=%s, TrackID=%d, Time=%.2f ps, Pos=(%.2f, %.2f, %.2f) Å\n",
                    //         name.Data(), trackID, time_ps, x_ang, y_ang, z_ang);
                    // }
                    g1->SetPoint(g1->GetN(), x_ang, z_ang);
                    numPointsPlotted++;
                    addedTrackIDsPerEvent[eventID].push_back(trackID);
                }
            }
            
            if (g1->GetN() > 0) {
                if (isDamaged) {
                    graphsDamaged.push_back(g1);
                } else {
                    graphsNormal.push_back(g1);
                }
            } else {
                delete g1;
            }
        }
        
        for (auto g : graphsNormal) g->Draw("P Same");
        for (auto g : graphsDamaged) g->Draw("P Same");
        
    }
    c->Update();
    
    std::string res;
    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        c->SaveAs(outFilename.Data());
    }
}
