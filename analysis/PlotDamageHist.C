/*
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <TString.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <vector>

void PlotDamageHist(TString filename) {
    if (filename == "") {
        printf("Please provide a valid ROOT filename.\n");
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
    TTree *tree = (TTree*)file->Get("Damage");
    if (!tree) {
        printf("Error: tree not found in file %s\n", filename.Data());
        return;
    }

    // Bind variables
    Int_t    eventID;
    Int_t    modelID;
    char     chainID[256];
    Int_t    residueID;
    char     compoundName[256];
    char     atomOrMoleculeName[256];
    Int_t    damageType;
    tree->SetBranchAddress("EventId", &eventID);
    tree->SetBranchAddress("ModelId", &modelID);
    tree->SetBranchAddress("ChainId", &chainID);
    tree->SetBranchAddress("ResidueId", &residueID);
    tree->SetBranchAddress("CompoundName", compoundName);
    tree->SetBranchAddress("AtomOrMoleculeName", atomOrMoleculeName);
    tree->SetBranchAddress("DamageType", &damageType);

    std::vector<TString> chainIDs;
    for (Long64_t i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        TString chainIDStr(chainID);
        if (std::find(chainIDs.begin(), chainIDs.end(), chainIDStr) == chainIDs.end()) {
            chainIDs.push_back(chainIDStr);
        }
    }
    std::cout << "Found " << chainIDs.size() << " unique chain IDs." << std::endl;
    if (chainIDs.size() == 0) {
        std::cout << "No chain IDs found in the Damage tree." << std::endl;
        return;
    }
    TString chainID1 = chainIDs[0];
    TString chainID2 = chainIDs.size() > 1 ? chainIDs[1] : "";

    Int_t nBasePair;
    std::cout << "Enter number of base pairs in DNA: ";
    std::cin >> nBasePair;

    TCanvas *c1 = new TCanvas("c1", "Base-Pair Damage Distribution", 600, 600);
    TH1I *hBpDamageDist = new TH1I("hBpDamageDist", ";Base-Pair ID;Number of Damages", nBasePair, 1, nBasePair+1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.07); 

    for (int i=0; i<tree->GetEntries(); i++) {
        if (residueID < 0 || residueID > nBasePair) std::cout << "hi" << std::endl; // Skip invalid residue IDs
        tree->GetEntry(i);
        Int_t basePairID = residueID;
        if (TString(chainID) == chainID2) 
            basePairID = nBasePair - residueID + 1; // Convert to base-pair ID
        hBpDamageDist->Fill(basePairID);
    }
    // Draw the histogram
    hBpDamageDist->Draw();
    std::string res;
    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        hBpDamageDist->SaveAs(outFilename);
    }

    TCanvas *c2 = new TCanvas("c2", "Nucleotide Damage Distribution", 600, 600);
    TH1I *hNucDamageDist = new TH1I("hNucDamageDist", ";Nucleotide ID;Number of Damages", nBasePair*2, 1, nBasePair*2+1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.07); 

    for (int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        Int_t nucleotideId = residueID;
        if (TString(chainID) == chainID2)
            nucleotideId += nBasePair;         // Shift nucleotide ID for second chain
        hNucDamageDist->Fill(nucleotideId);
    }
    hNucDamageDist->Draw();

    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        hNucDamageDist->SaveAs(outFilename.Data());
    }
}