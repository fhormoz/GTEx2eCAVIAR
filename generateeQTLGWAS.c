#include <cmath>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>

#include <stdlib.h>     /* atof */
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
	if(string(dirp->d_name).find(".peak")!=string::npos)
    		files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

void generateeQTLGWAS(map <string, vector <string> > & gtexCisSig, vector <string> snpIDs, vector <double> zScoreGWAS, string outputFileName, bool signOrNoSign) {	
	ofstream outputStream;
	
	for(int i = 0; i < snpIDs.size(); i++) {
		string snpId = snpIDs[i];
		vector <string> cisSig = gtexCisSig[snpId];
		for(int j = 0; j < cisSig.size(); j++) {
			string tmp = cisSig[j];
			string tmp2;
			string tmp3;
			string geneId = tmp.substr(0, tmp.find("\t"));
			outputStream.open( (outputFileName+ geneId + ".eqtl").c_str() , std::ofstream::out | std::ofstream::app);
			//TODO not sure which one is the z-score
			outputStream << snpId << "\t";
			tmp2 = tmp.substr(tmp.find("\t")+1);
			tmp3 = tmp2.substr(tmp2.find("\t")+1);
			double zscoretmp = atof(tmp3.substr(0, tmp3.find("\t")).c_str());
			if(signOrNoSign)	outputStream << zscoretmp << endl;  
			else			outputStream << abs(zscoretmp) << endl;	
			outputStream.close();
			outputStream.open( (outputFileName + geneId + ".gwas").c_str(), std::ofstream::out | std::ofstream::app);
			outputStream << snpId << "\t";
			if(signOrNoSign)	outputStream << zScoreGWAS[i] << endl;
			else 			outputStream << abs(zScoreGWAS[i]) << endl;
			outputStream.close();
		}	
	}	
}

void generateALLeQTLGWAS(map <string, vector <string> > & gtexCisSig, string dicName, string outputFileName, bool signOrNoSign) {
	vector <string> files;
        getDir(dicName, files);

	string line;
	ifstream test;
	vector <string> snpIDs;
	vector <double> zScoreGWAS;
	for(int i = 0; i < files.size(); i++) {
		cout << files[i] << endl;
		test.open((dicName + "/" + files[i]).c_str());
		while(std::getline(test, line)) {
			int indexFirstTab = line.find("\t");
			string snpID  = line.substr(0, indexFirstTab);
			double zScore = atof(line.substr(line.find("\t")+1, line.length()).c_str());
			snpIDs.push_back(snpID);
			zScoreGWAS.push_back(zScore);
		}
		string tmpfileName = files[i].substr(0, files[i].find(".peak"));
		generateeQTLGWAS(gtexCisSig, snpIDs, zScoreGWAS, outputFileName + "/" + tmpfileName, signOrNoSign);	
		snpIDs.erase( snpIDs.begin(), snpIDs.end() );
		zScoreGWAS.erase( zScoreGWAS.begin(), zScoreGWAS.end() );
		test.close();
	}
}
	
void obtainALLSNPsFromDirectory(string dicName, map <string, double> & hashMapSNP){										
	vector <string> files;
	getDir(dicName, files);

	string line;
	ifstream snpStream;
        for(int i = 0; i < files.size(); i++) {
		snpStream.open( (dicName + files[i]).c_str() );
		while(std::getline(snpStream, line)) {
			string snpID  = line.substr(0, line.find("\t"));
			float zScore = atof(line.substr(line.find("\t")+1, line.length()).c_str());
			hashMapSNP[snpID] = zScore;
		}
		snpStream.close();
	}	
}

void obtainGTExCisSignificant(string dosageFileName, map <string, vector <string> > & gtexCisSig, map <string, double> hashMapSNP2ZGWAS, vector <string> & geneGTExCisSig) {
	string line;
	map <string, int> geneMap;
	ifstream dosageStream;
	dosageStream.open(dosageFileName.c_str());
        
	while(std::getline(dosageStream, line)) {
		int indexFirstTab = line.find("\t");
                string snpID  = line.substr(0, indexFirstTab);
		string geneId = line.substr(indexFirstTab+1, line.find("\t", indexFirstTab+1)-indexFirstTab );
                if(hashMapSNP2ZGWAS.find(snpID) != hashMapSNP2ZGWAS.end()){
                        gtexCisSig[snpID].push_back(line.substr(line.find("\t")+1));
			if(geneMap.find(geneId) == geneMap.end()) {
				geneGTExCisSig.push_back(geneId);
				geneMap[geneId] = 0;
			}
		}
        }
}

int main(int argc, char *argv[]) {
	int oc = 0;
	int signValue = 0;
	bool signOrNoSign = true;
	string tissueName;
	string gwasName;
	string currentPath;
	map <string, double> hashMapSNP2ZGWAS;
	map <string, vector <string> > gtexCisSig;
	vector <string> geneGTExCisSig;

	vector <string> snpIDs;
	vector <double> zScoreGWAS;
	
	string cisGTEXFileName;
	
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
				cout << " -f GTEXCisFILE" << " specify the GTEx Cis significant folder" << endl;
				cout << " -s SIGN " << "specify if we need sign or not for Z-score and LD, 0 (everything is positive), 1 (sign is important)" << endl;
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
				cisGTEXFileName = string(optarg) + "/" + tissueName + "_Analysis.cis.eqtl"; 	
				break;
			case 's':
				signValue = atoi(optarg);
				signOrNoSign = ((signValue==0) ? false : true);
				break;
			case ':':
                        case '?':
                        default:
                                cout << "Strange" << endl;
                                break; 
		}
	}
	cout << cisGTEXFileName << endl;
	cout << currentPath + "/peak/" + gwasName + "/" + tissueName + "/" << endl;
	cout << currentPath + "/in/" + gwasName + "/" + tissueName << endl;
			
	obtainALLSNPsFromDirectory(currentPath + "/peak/" + gwasName + "/" + tissueName + "/" , hashMapSNP2ZGWAS);
	obtainGTExCisSignificant(cisGTEXFileName, gtexCisSig, hashMapSNP2ZGWAS, geneGTExCisSig);
	generateALLeQTLGWAS(gtexCisSig, currentPath + "/peak/" + gwasName + "/" + tissueName + "/" , currentPath + "/in/" + gwasName + "/" + tissueName + "/", signOrNoSign);	
	
}	
