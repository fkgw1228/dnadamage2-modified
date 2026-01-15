#include "DamageAnalyzer.hh"
#include <iostream>
#include <fstream>

void OutputResultsToFile(const std::vector<DamageAnalyzer>& analyzers, const std::string& outputFileName);

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file_name1> <input_file_name2> ... <output_file_name>" << std::endl;
        return 1;
    }

    std::vector<std::string> inputFileNames;
    std::string outputFileName = argv[argc-1];
    for (int i=1; i<argc-1; ++i) inputFileNames.push_back(argv[i]);

    std::vector<DamageAnalyzer> analyzers;
    try {
        for (std::string inputFileName : inputFileNames) {
            DamageAnalyzer analyzer(inputFileName);
            analyzer.Analyze();
            analyzers.push_back(analyzer);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during analysis: " << e.what() << std::endl;
        return 1;
    }
    // ask user confirmation before overwriting output file
    do {
        std::cout << "Output file '" << outputFileName << "' will be overwritten. Continue? (y/n): ";
        char response;
        std::cin >> response;
        if (response == 'n' || response == 'N') {
            std::cout << "Operation cancelled by user." << std::endl;
            return 0;
        } else if (response == 'y' || response == 'Y') {
            break;
        }
    } while (true);

    // Overwrite results to new output file
    OutputResultsToFile(analyzers, outputFileName);

    return 0;
}

void OutputResultsToFile(const std::vector<DamageAnalyzer>& analyzers, const std::string& outputFileName) {
    // Output file format:
    std::ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file " + outputFileName);
    }

    for (const auto& analyzer : analyzers) {
        outputFile << analyzer.GetLETMean() << " " << analyzer.GetLETStdDev() << " "
                   << analyzer.GetSSBToDSBRatio() << " " << 0 << "\n";
    }

    outputFile.close();
}