import os
import re
import time
import random
import glob
import numpy as npy
import sys, traceback
import subprocess
from myconfig import *
from subprocess import Popen, PIPE

def simulateData(zeQTLdata, zGWASdata, LD, causalSNP):
	_,snpCount = npy.shape(LD);
	true_casual = [0] * snpCount;
	true_casual[causalSNP] = zGWASdata[causalSNP];
	Zg = npy.random.multivariate_normal((npy.matrix(true_casual) * LD).tolist()[0], LD);
        true_casual[causalSNP] = zeQTLdata[causalSNP];
       	Ze = npy.random.multivariate_normal((npy.matrix(true_casual) * LD).tolist()[0], LD);
	return Ze, Zg;

def saveSimData(gtexSNPIds, simZeQTL, simZGWAS, outputfileZe, outputfileZg):
	fileOutZg = open(outputfileZg, "w");
	fileOutZe = open(outputfileZe, "w");
	for rsID, Zg in zip(gtexSNPIds, simZGWAS) :
		fileOutZg.write(rsID + "\t" + str(Zg)+"\n");
	for rsID, Ze in zip(gtexSNPIds, simZeQTL) :
                fileOutZe.write(rsID + "\t" + str(Ze)+"\n");
	fileOutZg.close();
	fileOutZe.close();

def maxSubmitReached(max):
	p1 = Popen(["qstat", "-u", "fhormoz"], stdout=PIPE)
	p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
	output = p2.communicate()[0]	
	if(int(output) < max) :
		return False;
	else:
		return True;

currentPath = os.getcwd();
gtex2rsId  = npy.genfromtxt(gwasFile, usecols=(3,4,), dtype='str');
print "Tissue", "rsID", "geneID", "CLPP" , "Adj-CLPP"; 
for tissue in glob.glob(currentPath+"/tmp/in/" + gwasName + "/*") :
	for eqtlFile in glob.glob(tissue + "/*.eqtl"):
		tissueName = tissue.replace(currentPath+"/tmp/in/" + gwasName + "/", '');
		eqtlName = eqtlFile.replace(currentPath + "/tmp/in/" + gwasName + "/" + tissueName + "/",'').replace('.eqtl','');	
		prefix_match = re.match(r"(.*_.*)_(rs\d*)_(.*)(ENSG.*)", eqtlName)
                if(prefix_match):
			rsID  = prefix_match.group(2),
			genID = prefix_match.group(4);	
		tmpGwasFile   = eqtlFile.replace('.eqtl', '.gwas'); 
		
		gtexSNPIds = npy.genfromtxt(eqtlFile, usecols=(0,), dtype='str'); 
		zeQTLdata  = npy.genfromtxt(eqtlFile, usecols=(1,));
		zGWASdata  = npy.genfromtxt(tmpGwasFile, usecols=(1,));
		indexrsId2gtexid = npy.where(rsID==gtex2rsId[:,1]);
		causalSNP = npy.where(gtexSNPIds==gtex2rsId[indexrsId2gtexid[0],0]);
		if(npy.ndim(zeQTLdata) == 0):
			continue;
		if( max(abs(zeQTLdata)) < 4 ):
			continue;
		if( max(abs(zGWASdata)) < 4 ):
                        continue;
		if (len(causalSNP[0]) ==0):
                        continue;
		#print currentPath + "/output/" + "/" + gwasName + "/sim/" + eqtlName +  "/*_col";	
		#/u/project/zarlab/fhormoz/code/GTEx_GWAS2/output/FG_100/sim/FG_100_rs560887_Breast_Mammary_TissueENSG00000152254.6/*_col
		realCOLVal = npy.genfromtxt(currentPath + "/output/" + "/" + gwasName + "/" + eqtlName + ".out_col", usecols=(1,));
		simCOLLists = [];
		for simFiles in glob.glob(currentPath + "/output/" + gwasName + "/sim/" + eqtlName +  "/*_col"):
			try:
				simCOLVal = npy.genfromtxt(simFiles, usecols=(1,));
			except ValueError:
				sys.stderr.write("fatal error" + simFiles);
				continue;
			except  IOError:
				sys.stderr.write("fatal error" + simFiles);
                                continue;	
			simCOLLists.append(simCOLVal[causalSNP]);
		countPasCut = sum(i<realCOLVal[causalSNP][0] for i in simCOLLists);
		if (len(simCOLLists) == 0):
			continue;
		print tissueName, rsID[0], genID,realCOLVal[causalSNP][0],float(countPasCut * 1.0)/float(len(simCOLLists));
