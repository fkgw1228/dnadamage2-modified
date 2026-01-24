#include "DamageAnalyzer.hh"
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

DamageAnalyzer::DamageAnalyzer(std::string inputFileName)
    : fInputFileName(inputFileName) {
        std::ifstream inputFile(fInputFileName);
        if (!inputFile.is_open()) {
            throw std::runtime_error("Could not open input file " + fInputFileName);
        }
}

void DamageAnalyzer::Analyze() {
    // Get first row showing LET value
    ParseLETValues(fInputFileName);

    // Parse damage records
    std::ifstream inputFile(fInputFileName);
    ParseDamageRecords(inputFile, fDamageRecords);
    inputFile.close();
    
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

void DamageAnalyzer::ParseLETValues(std::string inputFileName) {
    std::ifstream inputFile(inputFileName);
    std::string line;
    std::getline(inputFile, line);
    fLETMean = ParseLETMean(line);
    fLETStdDev = ParseLETStdDev(line);
    inputFile.close();

    return;
}

double DamageAnalyzer::ParseLETMean(const std::string& line) {
    // Example line format: "# LET = 10.5 +- 0.5 keV / um"
    size_t pos1 = line.find("LET = ");
    size_t pos2 = line.find(" +- ");
    if (pos1 == std::string::npos || pos2 == std::string::npos) {
        throw std::runtime_error("Invalid LET line format");
    }
    // Extract LET mean value
    std::string letValueStr = line.substr(pos1 + 6, pos2 - (pos1 + 6));
    return std::stod(letValueStr);
}

double DamageAnalyzer::ParseLETStdDev(const std::string& line) {
    // Example line format: "# LET = 10.5 +- 0.5 keV / um"
    size_t pos1 = line.find("+- ");
    size_t pos2 = line.find(" keV");
    if (pos1 == std::string::npos || pos2 == std::string::npos) {
        throw std::runtime_error("Invalid LET line format");
    }
    // Extract LET standard deviation value
    std::string letStdDevStr = line.substr(pos1 + 3, pos2 - (pos1 + 3));
    return std::stod(letStdDevStr);
}

void DamageAnalyzer::ParseDamageRecords(std::ifstream& inputFile, std::vector<DamageRecord>& records) {
    std::string line;
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty lines or comments
        std::istringstream iss(line);
        DamageRecord record;
        iss >> record.eventID >> record.dnaID >> record.chainID >> record.nucleotideId 
            >> record.compoundName >> record.atomOrMoleculeName >> record.damageType;
        records.push_back(record);
    }
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
    // Group damages by eventID, chainID, nucleotideId, compoundName
    // If there is at least one damage of the specified type, count as SSB
    for (const auto& record : records) {
        SSBRecord ssb{record.eventID, record.chainID, record.nucleotideId, record.damageType};
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
                    abs(ssb1.nucleotideId - ssb2.nucleotideId) <= bpDistance) {
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
    std::cout << std::setw(35) << "SSB to DSB Ratio(%): "                 << "   " << fSSBToDSBRatio*100        << std::endl;
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