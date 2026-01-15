#include "DamageAnalyzer.hh"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_damage_file>" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    DamageAnalyzer analyzer(inputFileName);

    try {
        analyzer.Analyze();
        analyzer.PrintResults();
    } catch (const std::exception& e) {
        std::cerr << "Error during analysis: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}