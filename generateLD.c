#include <cmath>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>

#include <stdlib.h>
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

double computeAvg(vector <double> data) {
	double total = 0;
	for(int i = 0; i < data.size(); i++)
		total = total + data[i];
	return(total/data.size());
}

double computeVar(vector <double> data, double mean) {
	double total = 0;
	for(int i = 0; i < data.size(); i++)
		total = total + ((data[i]-mean) * (data[i]-mean));
	return( total/data.size() );
}

double computeCov(vector <double> data1, vector <double> data2) {
	double mean1 = computeAvg(data1);
	double mean2 = computeAvg(data2);
	double total = 0;
	for(int i = 0; i < data1.size(); i++)
                total = total + ((data1[i]-mean1) * (data2[i]-mean2));
	return(total / data1.size());
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
	double mean1 = computeAvg(dataV1); 
	double var1  = computeVar(dataV1, mean1);	
	double mean2 = computeAvg(dataV2);
	double var2  = computeVar(dataV2, mean2);
	double cov   = computeCov(dataV1, dataV2);
	return (cov/(sqrt(var1)*sqrt(var2)));
}

void generateLD(map <string, string> & snpDosageData, vector <string> snpIDs, string outputFileName, bool signOrNoSign) {	
	vector < vector<double> > LD;
	ofstream outputStream;

	LD.resize(snpIDs.size());
	outputStream.open(outputFileName.c_str());	
	for(int i = 0; i < snpIDs.size(); i++){
		LD[i].resize(snpIDs.size());
		for(int j = 0; j < snpIDs.size(); j++) {
			LD[i][j] = computeCor(snpDosageData[snpIDs[i]], snpDosageData[snpIDs[j]]);
			if(signOrNoSign)	outputStream << LD[i][j] << " ";
			else 			outputStream << abs(LD[i][j]) << " ";
		}
		outputStream << endl;
	}
	outputStream.close();	
}

void generteALLLD(string dicName, map <string, string> & snpDosageData, bool signOrNoSign) {
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
		string tmpName = eqtlFiles[i].substr(0, eqtlFiles[i].find(".eqtl"));
		generateLD(snpDosageData, tmpSNPs, dicName+ tmpName + ".LD" , signOrNoSign);
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
	int signValue = 0;
	bool signOrNoSign = true;
	string tissueName;
	string gwasName;
	string currentPath;
	string genotypeFile;
	map <string, int> hashMapSNP;
	vector <string> snpLists;
	map <string, string> snpDosageData;

	 while ((oc = getopt(argc, argv, "vhl:f:o:t:g:s:")) != -1) {
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
	                        cout << " -s SIGN" << "specify the sign of LD, 0 (No sign, eveything is positive), 1 (Sing is important)" << endl;     
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
			case 's':
				signValue = atoi(optarg);
				signOrNoSign = ((signValue==1)? true : false);
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
	generteALLLD(currentPath + "/in/" + gwasName + "/" + tissueName + "/", snpDosageData, signOrNoSign);	
}
