#include <cmath>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

using namespace std;

int getDir (string dir, vector<string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {
		if (string(dirp->d_name).find("eqtl") != string::npos)
                	files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

double computeCor(string data1, string data2) {
	double d;
	vector <double> dataV1;
	vector <double> dataV2;
	stringstream ssin1(data1);
	stringstream ssin2(data2);
	while(!ssin1.eof()) {	
		ssin1 >> d;
		dataV1.push_back(d);
		ssin2 >> d;
		dataV2.push_back(d);
	}
	double sumProduct = 0;
	double sumV1Pow2  = 0;
	double sumV2Pow2  = 0;	
	for(int i = 0; i < dataV1.size(); i++) {
		sumProduct += dataV1[i] * dataV2[i];
		sumV1Pow2  += (dataV1[i]*dataV1[i]);
		sumV2Pow2  += (dataV2[i]*dataV2[i]);
	}	
	return( sumProduct/(sqrt(sumV1Pow2)*sqrt(sumV2Pow2)) );	
}

void generateLD(map <string, string> & snpDosageData, vector <string> snpIDs, string outputFileName) {	
	vector < vector<double> > LD;
	ofstream outputStream;

	LD.resize(snpIDs.size());
	outputStream.open(outputFileName.c_str());	
	for(int i = 0; i < snpIDs.size(); i++){
		LD[i].resize(snpIDs.size());
		for(int j = 0; j < snpIDs.size(); j++) {
			LD[i][j] = computeCor(snpDosageData[snpIDs[i]], snpDosageData[snpIDs[j]]);
			outputStream << LD[i][j] << " ";
		}
		outputStream << endl;
	}
	outputStream.close();	
}

void generteALLLD(string dicName, map <string, string> & snpDosageData) {
	vector <string> eqtlFiles;
	vector <string> tmpSNPs;
        getDir(dicName, eqtlFiles);
		
	string line;                  
        ifstream snpStream;           
        for(int i = 0; i < eqtlFiles.size(); i++) {
                snpStream.open( (dicName + eqtlFiles[i]).c_str() );
                while(std::getline(snpStream, line))    {
                        string snpID = line.substr(0, line.find("\t"));
                        tmpSNPs.push_back(snpID);
                }
		string tmpName = eqtlFiles[i].substr(0, eqtlFiles[i].find_last_of(".eqtl"));
		generateLD(snpDosageData, tmpSNPs, dicName+ tmpName + ".LD" );
                snpStream.close();
		tmpSNPs.erase(tmpSNPs.begin(), tmpSNPs.end());
        }
}
	
void obtainALLSNPsFromDirectory(string dicName, map <string, int> & hashMapSNP){										
	vector <string> eqtlFiles;
	getDir(dicName, eqtlFiles);

	string line;
	ifstream snpStream;
        for(int i = 0; i < eqtlFiles.size(); i++) {
		snpStream.open( (dicName + eqtlFiles[i]).c_str() );
		while(std::getline(snpStream, line))    {
			string snpID = line.substr(0, line.find("\t"));
			hashMapSNP[snpID] = 0;
		}
		snpStream.close();
	}	
}

void obtainSNPsFromFile(string snpFileName, map <string, int> & hashMapSNP, vector <string> & snpLists) {
	int lineNum = 0;
	string line;
	ifstream snpStream;
        snpStream.open(snpFileName.c_str());
	while(std::getline(snpStream, line))    {
                string snpID = line.substr(0, line.find("\t"));
		snpLists.push_back(snpID);
                hashMapSNP[snpID] = lineNum;
                lineNum++;
        }	
}

void obtainDosageData(string dosageFileName, map <string, string> & snpDosageData, map <string, int> hashMapSNP) {
	string line;
	ifstream dosageStream;
	dosageStream.open(dosageFileName.c_str());
	//remove the GTEx header
	std::getline(dosageStream, line);     
        while(std::getline(dosageStream, line))    {
                string snpID = line.substr(0, line.find("\t"));
                if(hashMapSNP.find(snpID) != hashMapSNP.end()){
                        snpDosageData[snpID] = line.substr(line.find("\t")+1);
                }
        }
}

int main( int argc, char *argv[] ) {
	int oc = 0;
	string tissueName;
	string gwasName;
	string currentPath;
	string genotypeFile;
	map <string, int> hashMapSNP;
	vector <string> snpLists;
	map <string, string> snpDosageData;

	 while ((oc = getopt(argc, argv, "vhl:f:o:t:g:")) != -1) {
                switch (oc) {
                        case 'v':
                                cout << "version 0.0:" << endl;
                        case 'h':
                                cout << "Options: " << endl;
                                cout << " -h help"   << endl;
                                cout << " -t TISSUENAME"  << " specify the name of tissue" << endl;
                                cout << " -g GWASNAME"    << " specify the name of GWAS"   << endl;
                                cout << " -o TMPFILE"     << " specify the tmp folder"     << endl;
                                cout << " -f Genotype File" << " specify the GTEx Genotype folder" << endl;
                                return(0);
                        case 't':
                                tissueName  = string(optarg);
                                break;
                        case 'g':
                                gwasName    = string(optarg);
                                break;
                        case 'o':
                                currentPath = string(optarg);
                                break;
                        case 'f':
                                genotypeFile = string(optarg) + "/" + tissueName + "_Analysis.snps.txt";
                                break;
                        case ':':
                        case '?':
                        default:
                                cout << "Strange" << endl;
                                break;
                }
        }

	obtainALLSNPsFromDirectory(currentPath + "/in/" + gwasName + "/" + tissueName + "/", hashMapSNP);
	obtainDosageData(genotypeFile, snpDosageData, hashMapSNP);
	generteALLLD(currentPath + "/in/" + gwasName + "/" + tissueName + "/", snpDosageData);	
}
