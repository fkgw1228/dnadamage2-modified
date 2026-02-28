/*
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/


#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

struct DamageRecord {
    int eventID;
    int modelID;
    std::string chainID;
    int residueId;
    std::string compoundName;
    std::string atomOrMoleculeName;
    int damageType;         // 1 for Direct Damage, 2 for Indirect Damage
};

struct SSBRecord {
    int eventID;
    std::string chainID;
    int residueId;
    int damageType;         // 1 for Direct Damage, 2 for Indirect Damage
    bool operator==(const SSBRecord& other) const {
        return eventID == other.eventID
            && chainID == other.chainID
            && residueId == other.residueId;
    }
};

class DamageAnalyzer {
public:
    // Constructor
    DamageAnalyzer(std::string inputFileName);

    // Main analysis function
    void   Analyze();
    void   PrintResults();
    
    // Getters
    std::string GetInputFileName() const { return fInputFileName; }
    double GetLETMean() const { return fLETMean; }
    double GetLETStdDev() const { return fLETStdDev; }
    int    GetTotalDirectDamages() const { return fTotalDirectDamages; }
    int    GetTotalIndirectDamages() const { return fTotalIndirectDamages; }
    int    GetDirectSSBs() const { return fDirectSSBs; }
    int    GetIndirectSSBs() const { return fIndirectSSBs; }
    int    GetTotalSSBs() const { return fTotalSSBs; }
    int    GetTotalDSBs() const { return fTotalDSBs; }
    double GetIndirectToTotalRatio() const { return fIndirectToTotalRatio; }
    double GetSSBToDSBRatio() const { return fSSBToDSBRatio; }

private:
    void   ParseLETValues();
    void   ParseDamageRecords();

    void   NoResult();

    int    CountDamages(const std::vector<DamageRecord>& records, int damageType);
    int    CountSSBs(const std::vector<DamageRecord>& records, std::vector<SSBRecord>& ssbRecords);
    int    CountDSBs(const std::vector<SSBRecord>& ssbRecords, int bpDistance = 10);
    int    CountDirectSSBs(const std::vector<SSBRecord>& ssbRecords);
    int    CountIndirectSSBs(const std::vector<SSBRecord>& damageRecords);

    int fBasePairNum = 0;
    std::string fInputFileName;
    double fLETMean, fLETStdDev;
    int fTotalDirectDamages, fTotalIndirectDamages, fTotalDamages;
    int fDirectSSBs, fIndirectSSBs, fTotalSSBs, fTotalDSBs;
    double fIndirectToTotalRatio, fSSBToDSBRatio;
    std::vector<DamageRecord> fDamageRecords;
    std::vector<SSBRecord> fSSBRecords;
};

DamageAnalyzer::DamageAnalyzer(std::string inputFileName)
    : fInputFileName(inputFileName) {
        // Check if file exists
        TFile* f = TFile::Open(fInputFileName.c_str());
        if (!f || f->IsZombie()) {
             throw std::runtime_error("Could not open input file " + fInputFileName);
        }
        f->Close();

        std::cout << "Enter the base pair number of the DNA to be analyzed: ";
        std::cin >> fBasePairNum;
}

void DamageAnalyzer::Analyze() {
    // Get LET values
    ParseLETValues();

    // Parse damage records
    ParseDamageRecords();
    
    if (fDamageRecords.empty()) {
        NoResult();
        return;
    }
    // Analyze damage records
    fTotalDirectDamages =   CountDamages(fDamageRecords, 1);
    fTotalIndirectDamages = CountDamages(fDamageRecords, 2);
    fTotalDamages = fTotalDirectDamages + fTotalIndirectDamages;
    fIndirectToTotalRatio = fTotalDamages > 0 ? fTotalIndirectDamages / (double)(fTotalDamages) : 0.0;
    
    // SSB calculation
    fTotalSSBs = CountSSBs(fDamageRecords, fSSBRecords);
    // Direct and Indirect SSBs
    fDirectSSBs = CountDirectSSBs(fSSBRecords);
    fIndirectSSBs = CountIndirectSSBs(fSSBRecords);
    // DSB calculation
    fTotalDSBs = CountDSBs(fSSBRecords);
    fSSBToDSBRatio = (fTotalDSBs > 0) ? (fTotalSSBs / (double)(fTotalDSBs)) : 0.0;

    return;
}

void DamageAnalyzer::ParseLETValues() {
    TFile* f = TFile::Open(fInputFileName.c_str());
    if (!f || f->IsZombie()) {
        throw std::runtime_error("Could not open ROOT file " + fInputFileName);
    }
    TTree* t = (TTree*)f->Get("LET");
    if (!t) {
        // Maybe default file without LET info?
        // Or handle error
        fLETMean = 0.0;
        fLETStdDev = 0.0;
        f->Close();
        return;
    }
    Double_t mean, std;
    t->SetBranchAddress("meanLET", &mean);
    t->SetBranchAddress("stdLET", &std);
    t->GetEntry(0);
    fLETMean = mean;
    fLETStdDev = std;
    f->Close();
}

void DamageAnalyzer::ParseDamageRecords() {
    TFile* f = TFile::Open(fInputFileName.c_str());
    if (!f || f->IsZombie()) {
        throw std::runtime_error("Could not open ROOT file " + fInputFileName);
    }
    TTree* t = (TTree*)f->Get("Damage");
    if (!t) {
        f->Close();
        return;
    }
    
    Int_t eventId, modelId, residueId, damageType;
    char chainId[128], compoundName[128], atomName[128];
    
    t->SetBranchAddress("EventId", &eventId);
    t->SetBranchAddress("ModelId", &modelId);
    t->SetBranchAddress("ChainId", chainId);
    t->SetBranchAddress("ResidueId", &residueId);
    t->SetBranchAddress("CompoundName", compoundName);
    t->SetBranchAddress("AtomOrMoleculeName", atomName);
    t->SetBranchAddress("DamageType", &damageType);
    
    Long64_t nEntries = t->GetEntries();
    for(Long64_t i=0; i<nEntries; ++i) {
        t->GetEntry(i);
        DamageRecord record;
        record.eventID = eventId;
        record.modelID = modelId;
        record.chainID = std::string(chainId);
        record.residueId = residueId;
        record.compoundName = std::string(compoundName);
        record.atomOrMoleculeName = std::string(atomName);
        record.damageType = damageType;
        fDamageRecords.push_back(record);
    }
    f->Close();
}

int DamageAnalyzer::CountDamages(const std::vector<DamageRecord>& records, int damageType) {
    int count = 0;
    for (const auto& record : records) {
        if (record.damageType == damageType) {
            count++;
        }
    }
    return count;
}

int DamageAnalyzer::CountSSBs(const std::vector<DamageRecord>& records, std::vector<SSBRecord>& ssbRecords) {
    // Group damages by eventID, chainID, residueId, compoundName
    // If there is at least one damage of the specified type, count as SSB
    for (const auto& record : records) {
        SSBRecord ssb{record.eventID, record.chainID, record.residueId, record.damageType};
        if (std::find(ssbRecords.begin(), ssbRecords.end(), ssb) == ssbRecords.end()) {
            ssbRecords.push_back(ssb);
        }
    }
    int ssbCount = ssbRecords.size();

    return ssbCount;
}

int DamageAnalyzer::CountDSBs(const std::vector<SSBRecord>& ssbRecords, int bpDistance) {
    int dsbCount = 0;
    // If two SSBs are on opposite strands within bpDistance at the same event, count as DSB
    // Don't double count SSBs that already form DSBs

    // Group SSBs by eventID
    std::map<int, std::vector<SSBRecord>> eventSSBMap;
    for (const auto& ssb : ssbRecords) 
        eventSSBMap[ssb.eventID].push_back(ssb);

    for (const auto& eventAndSSBs : eventSSBMap) {
        int eventID = eventAndSSBs.first;
        const auto& ssbs = eventAndSSBs.second;
        // flags for SSB used in DSB
        std::vector<bool> used(ssbs.size(), false);

        for (size_t i=0; i<ssbs.size(); i++) {
            if (used[i]) continue;
            for (size_t j=i+1; j<ssbs.size(); j++) {
                if (used[j]) continue;
                
                const auto& ssb1 = ssbs[i];
                const auto& ssb2 = ssbs[j];
                if (ssb1.chainID != ssb2.chainID && 
                    abs(ssb1.residueId - abs(ssb2.residueId - fBasePairNum)) <= bpDistance) {
                    dsbCount++;
                    used[i] = true;
                    used[j] = true;
                    break; // Move to next SSB after forming a DSB
                }
            }
        }
    }

    return dsbCount;
}

int DamageAnalyzer::CountDirectSSBs(const std::vector<SSBRecord>& ssbRecords) {
    int directSSBCount = 0;
    for (const auto& ssb : ssbRecords) {
        if (ssb.damageType == 1) { // Direct Damage
            directSSBCount++;
        }
    }
    return directSSBCount;
}

int DamageAnalyzer::CountIndirectSSBs(const std::vector<SSBRecord>& damageRecords) {
    int indirectSSBCount = 0;
    for (const auto& record : damageRecords) {
        if (record.damageType == 2) { // Indirect Damage
            indirectSSBCount++;
        }
    }
    return indirectSSBCount;
}

void DamageAnalyzer::PrintResults() {
    std::cout << std::setw(35) << "File Name: "                           << "   " << fInputFileName            << std::endl;
    std::cout << std::setw(35) << "LET: "                                 << "   " << fLETMean << " +- " << fLETStdDev << " keV/um" << std::endl;
    std::cout << std::setw(35) << "Total Direct Damages: "                << "   " << fTotalDirectDamages       << std::endl;
    std::cout << std::setw(35) << "Total Indirect Damages: "              << "   " << fTotalIndirectDamages     << std::endl;
    std::cout << std::setw(35) << "Total SSBs: "                          << "   " << fTotalSSBs                << std::endl;
    std::cout << std::setw(35) << "Direct SSBs: "                         << "   " << fDirectSSBs               << std::endl;
    std::cout << std::setw(35) << "Indirect SSBs: "                       << "   " << fIndirectSSBs             << std::endl;
    std::cout << std::setw(35) << "Total DSBs: "                          << "   " << fTotalDSBs                << std::endl;
    std::cout << std::setw(35) << "Indirect to Total Damage Ratio:(%) "   << "   " << fIndirectToTotalRatio*100 << std::endl;
    std::cout << std::setw(35) << "SSB to DSB Ratio: "                    << "   " << fSSBToDSBRatio            << std::endl;
}

void DamageAnalyzer::NoResult() {
    fTotalDirectDamages = 0;
    fTotalIndirectDamages = 0;
    fDirectSSBs = 0;
    fIndirectSSBs = 0;
    fTotalSSBs = 0;
    fTotalDSBs = 0;
    fIndirectToTotalRatio = 0.0;
    fSSBToDSBRatio = 0.0;
    return;
}

void OutputResultsToFile(const std::vector<DamageAnalyzer>& analyzers, const std::string& outputFileName) {
    // Output file format:
    std::ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file " + outputFileName);
    }

    for (const auto& analyzer : analyzers) {
        outputFile << analyzer.GetLETMean() << " " << analyzer.GetSSBToDSBRatio() << " "
                   << analyzer.GetLETStdDev() << " " << 0 << "\n";
    }

    outputFile.close();
}

void PrintDamageSummary(std::vector<TString> inputFileNames) {
    if (inputFileNames.empty()) {
        std::cout << "Please specify the input damage file (.root)." << std::endl;
        std::cout << 
            "Usage example: .x PrintDamagesSummary({\"data/strand_breaks/StrandBreakInfo.txt\"})"
                << std::endl;
        return;
    }

    std::string outputFileName;
    std::string res;
    std::cout << "Do you want to output SSB/DSB ratio for each LET to a file? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        std::cout << "Enter the name of output file: ";
        std::cin >> outputFileName;
    }

    std::vector<DamageAnalyzer> analyzers;
    for ( TString inputFileName : inputFileNames) {
        DamageAnalyzer analyzer(inputFileName.Data());
        analyzer.Analyze();
        analyzer.PrintResults();
        analyzers.push_back(analyzer);
    }
    if (!outputFileName.empty()) {
        OutputResultsToFile(analyzers, outputFileName);
    }
}
