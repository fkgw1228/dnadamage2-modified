#include <string>
#include <vector>

struct DamageRecord {
    int eventID;
    int dnaID;
    std::string chainID;
    int nucleotideId;
    std::string compoundName;
    std::string atomOrMoleculeName;
    int damageType;         // 1 for Direct Damage, 2 for Indirect Damage
};

struct SSBRecord {
    int eventID;
    std::string chainID;
    int nucleotideId;
    int damageType;         // 1 for Direct Damage, 2 for Indirect Damage
    bool operator==(const SSBRecord& other) const {
        return eventID == other.eventID
            && chainID == other.chainID
            && nucleotideId == other.nucleotideId;
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
    void   ParseLETValues(std::string inputFileName);
    double ParseLETMean(const std::string& line);
    double ParseLETStdDev(const std::string& line);
    void   ParseDamageRecords(std::ifstream& inputFile, std::vector<DamageRecord>& records);
    int    CountDamages(const std::vector<DamageRecord>& records, int damageType);
    int    CountSSBs(const std::vector<DamageRecord>& records, std::vector<SSBRecord>& ssbRecords);
    int    CountDSBs(const std::vector<SSBRecord>& ssbRecords, int bpDistance = 10);
    int    CountIndirectSSBs(const std::vector<SSBRecord>& damageRecords);
    int    CountDirectSSBs(const std::vector<SSBRecord>& ssbRecords);
    void   NoResult();

    std::string fInputFileName;
    double fLETMean, fLETStdDev;
    int fTotalDirectDamages, fTotalIndirectDamages, fTotalDamages;
    int fDirectSSBs, fIndirectSSBs, fTotalSSBs, fTotalDSBs;
    double fIndirectToTotalRatio, fSSBToDSBRatio;
    std::vector<DamageRecord> fDamageRecords;
    std::vector<SSBRecord> fSSBRecords;
};