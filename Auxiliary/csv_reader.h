#ifndef __CSV_READER_H__
#define __CSV_READER_H__

#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

class csv_reader {
    
    
    
    
    // based on [https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/]
    
    // example usage:
    //    int main() {
    //        // Read three_cols.csv and ones.csv
    //        std::vector<std::pair<std::string, std::vector<double>>> three_cols = read_csv("three_cols.csv");
    //        std::vector<std::pair<std::string, std::vector<double>>> ones = read_csv("ones.csv");
    //
    //        return 0;
    //    }
    
public:
    csv_reader(){};
    ~csv_reader(){};
    
    
    std::vector<std::pair<std::string, double>> read_csv(std::string filename){
        // Reads a CSV file into a vector of <string, vector<int>> pairs where
        // each pair represents <column name, column values>
        
        // Create a vector of <string, int vector> pairs to store the result
        std::vector<std::pair<std::string, double>> result;
        
        // Create an input filestream
        std::ifstream myFile(filename);
        
        // Make sure the file is open
        if(!myFile.is_open()) throw std::runtime_error("Could not open file");
        
        // Helper vars
        std::string line, colname;
        std::string token;
        
        // read header line
        std::getline(myFile, line);
        
        // Read data, line by line
        while(std::getline(myFile, line)) {
            // Create a stringstream of the current line
            std::stringstream ss(line);
            // std::cout << ss.str() << std::endl;
            
            std::pair<std::string, double> cut;
            std::getline(ss, token, ',');
            cut.first = token;
            if (strcmp("",token.c_str())==0) {
                break;
            }
            std::getline(ss, token, ',');
            cut.second = (double)(std::stof( token ));
            
            
            result.push_back(cut);
        }
        
        // Close file
        myFile.close();
        
        return result;
    }
    
    
    
protected:
    
};
#endif

