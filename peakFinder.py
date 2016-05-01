import os
import optparse
import numpy as npy
import subprocess
import scipy.stats as st

def calculateP2Zscore(betaVec, betaSDVec, pvalVec):
	resVec = [0] * len(betaVec);
	index = 0;
	print betaVec;
	for beta, betasd, pval in zip(betaVec, betaSDVec, pvalVec):
		stat = abs(st.norm.ppf(pval/2));
		if(beta * betasd < 0):
			stat = -stat;
		resVec[index] = stat;
		index = index + 1;
	return(resVec);

def main(parser):
        (options, args) = parser.parse_args();
	tissueName = options.tissue;
	gwasName   = options.gwasname;	
	sigSNPName = options.sigsnpfile;
	gwasFile   = options.gwasfile;
	outfolder  = options.outfolder;
	gtexPath   = options.gtexfolder;
	
	chrom    = npy.loadtxt(open(gwasFile, 'rb'), usecols=(0,), dtype='str', skiprows=0);
	startPos = npy.loadtxt(open(gwasFile, 'rb'), usecols=(1,), dtype='int', skiprows=0);
	endPos   = npy.loadtxt(open(gwasFile, 'rb'), usecols=(2,), dtype='int', skiprows=0);
	gtexID   = npy.loadtxt(open(gwasFile, 'rb'), usecols=(3,), dtype='str', skiprows=0);
	rsID     = npy.loadtxt(open(gwasFile, 'rb'), usecols=(4,), dtype='str', skiprows=0);
	beta     = npy.loadtxt(open(gwasFile, 'rb'), usecols=(5,), dtype='float', skiprows=0);
	betaSD   = npy.loadtxt(open(gwasFile, 'rb'), usecols=(6,), dtype='float', skiprows=0);
	pVal     = npy.loadtxt(open(gwasFile, 'rb'), usecols=(7,), dtype='float', skiprows=0);
	#zScore	 = npy.loadtxt(open(gwasFile, 'rb'), usecols=(8,), dtype='float', skiprows=0);
	
	significantSNPsRSID = npy.loadtxt(open(sigSNPName, 'rb'), usecols=(0,), dtype='str', skiprows=0);
	newList = [];
	for snp in significantSNPsRSID:
		newList.append(npy.where(rsID==snp));

	#Generate ALL the needed folders
	if not os.path.exists(outfolder):
		os.makedirs(outfolder);
	if not os.path.exists(outfolder + "/peak/"):
                 os.makedirs(outfolder + "/peak/");
	if not os.path.exists(outfolder + "/in/"):
		 os.makedirs(outfolder + "/in/");
	if not os.path.exists(outfolder + "/peak/" + gwasName):
        	os.makedirs(outfolder + "/peak/"+ gwasName);
	if not os.path.exists(outfolder + "/peak/" + gwasName + "/" + tissueName):
		os.makedirs(outfolder + "/peak/" + gwasName + "/" + tissueName);
	if not os.path.exists(outfolder + "/in/" + gwasName):
                os.makedirs(outfolder + "/in/"+ gwasName);
        if not os.path.exists(outfolder + "/in/" + gwasName + "/" + tissueName):
                os.makedirs(outfolder + "/in/" + gwasName + "/" + tissueName);	
	
	count = 0;	
	for indexs in newList:
		peakOutFile = open(outfolder + "/peak/"+ gwasName  +"/"+ tissueName + "/" +gwasName+"_"+str(rsID[indexs][0])+"_"+ tissueName +".peak", 'w');
		rangeIndexStart = newList[count][0]-50;
		rangeIndexEnd   = newList[count][0]+50; 
		gtexIDAroundPeak = gtexID[rangeIndexStart:rangeIndexEnd];
		zScore = calculateP2Zscore(beta[rangeIndexStart:rangeIndexEnd], betaSD[rangeIndexStart:rangeIndexEnd], pVal[rangeIndexStart:rangeIndexEnd]);	
		for snpId, zVal in zip(gtexIDAroundPeak, zScore) :
			peakOutFile.write(snpId + "\t" + str(zVal) + "\n");
		peakOutFile.close();
		count = count + 1;
	print "We Generate the Peaks files";
	generateeQTLGWASScript = "./generateeQTLGWAS -t " + tissueName + " -g " + gwasName + " -o " + \
				 outfolder + " -f " + gtexPath + "/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/";
	print generateeQTLGWASScript;	
	subprocess.Popen(generateeQTLGWASScript ,shell=True).wait();
	print "We Generate the eQTL and GWAS files";
	generateLDScript       = "./generateLD -t " + tissueName + " -g " + gwasName + " -o " + \
				 outfolder + " -f " + gtexPath + "/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices/" ;
	print generateLDScript;
	subprocess.Popen(generateLDScript ,shell=True).wait();	
	
if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-n", "--gname", dest="gwasname",
                default="", type="string",
                help="specify the name of GWAS Study");
	parser.add_option("-g", "--gtexpath", dest="gtexfolder",
                default="", type="string",
                help="specify the Path to GTEx Folder");
        parser.add_option("-f", "--input", dest="gwasfile",
                default="", type="string",
                help="specify the full path to GWAS File")
	parser.add_option("-t", "--tissueName", dest="tissue", default="",
                type="string", help="the tissue name")
	parser.add_option("-s", "--sigSNPName", dest="sigsnpfile", default="",
                type="string", help="the file whcih have the significnat SNP rsID name")	
	parser.add_option("-o", "--outfolder", dest="outfolder", default="",
                type="string", help="the full path to generate files needed for eCAVIAR")

        main(parser);
