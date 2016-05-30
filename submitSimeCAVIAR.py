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
        if(zGWASdata[causalSNP]>0):
                true_casual[causalSNP] = max(zGWASdata[causalSNP],5);
        else:
                true_casual[causalSNP] = min(zGWASdata[causalSNP],-5);
        Zg = npy.random.multivariate_normal((npy.matrix(true_casual) * LD).tolist()[0], LD);
        if(zeQTLdata[causalSNP]>0):
                true_casual[causalSNP] = max(zeQTLdata[causalSNP], 5);
        else:
                true_casual[causalSNP] = min(zeQTLdata[causalSNP], -5);
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


TEMPLATE_SERIAL = """
#####################################
#$ -S /bin/bash
#$ -cwd
#$ -N {name}
#$ -e {errfile}
#$ -o {logfile}
#$ -pe make {slots}
#$ -l h_data=1G
#$ -l h_rt=15:00:00
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
cd {ecaviarPath}
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""

currentPath = os.getcwd();
outputPath  = currentPath + "/output/" + gwasName + "/sim/";

if not os.path.exists(currentPath):      
        os.makedirs(currentPath);
if not os.path.exists(currentPath+"/output/"):
        os.makedirs(currentPath+"/output/");
if not os.path.exists(currentPath+"/output/"+gwasName):
        os.makedirs(currentPath+"/output/"+gwasName);
if not os.path.exists(outputPath):
	os.makedirs(outputPath);
if not os.path.exists(outputPath+"/log/"):      
        os.makedirs(outputPath+"/log/");

gtex2rsId  = npy.genfromtxt(gwasFile, usecols=(3,4,), dtype='str');
for tissue in glob.glob(currentPath+"/tmp/in/" + gwasName + "/*") :
	for eqtlFile in glob.glob(tissue + "/*.eqtl"):
		tissueName = tissue.replace(currentPath+"/tmp/in/" + gwasName + "/", '');
		eqtlName = eqtlFile.replace(currentPath + "/tmp/in/" + gwasName + "/" + tissueName + "/",'').replace('.eqtl','');	
		prefix_match = re.match(r"(.*_.*)_(rs\d*)_(.*)(ENSG.*)", eqtlName)
                if(prefix_match):
			rsID  = prefix_match.group(2),
			genID = prefix_match.group(4);	
		ldFile        = eqtlFile.replace('.eqtl', '.LD');
		tmpGwasFile   = eqtlFile.replace('.eqtl', '.gwas'); 
		gtexSNPIds = npy.genfromtxt(eqtlFile, usecols=(0,), dtype='str'); 
		zeQTLdata  = npy.genfromtxt(eqtlFile, usecols=(1,));
		zGWASdata  = npy.genfromtxt(tmpGwasFile, usecols=(1,));
		indexrsId2gtexid = npy.where(rsID==gtex2rsId[:,1]);
		causalSNP = npy.where(gtexSNPIds==gtex2rsId[indexrsId2gtexid[0],0]);
		if(npy.ndim(zeQTLdata) == 0):
                        continue;

		LD = npy.loadtxt(ldFile);
                rowLD, colLD = LD.shape;
                LDdet = npy.linalg.det(LD);
                while (LDdet <= 0):
                        LD = LD + 0.1 * npy.eye(rowLD);
                        LDdet = npy.linalg.det(LD);
		if( max(abs(zeQTLdata)) < 4 ):
			continue;
		if( max(abs(zGWASdata)) < 4 ):
                        continue;
		if (len(causalSNP[0]) ==0):
			continue;		
		if not os.path.exists(outputPath + "/" + eqtlName):
                        os.makedirs(outputPath + "/" + eqtlName);	
		script = "";
		errfile      = outputPath + "/log/Sim" + eqtlName  + ".err"
		logfile      = outputPath + "/log/Sim" + eqtlName  + ".log"
		scriptfile   = outputPath + "/log/Sim" + eqtlName;
		for index in range(400):
			#we make the index of rsID SNP to be causal.
			simZeQTL, simZGWAS = simulateData(zeQTLdata, zGWASdata, LD, causalSNP[0]);
	
			outputfile   = outputPath + "/" + eqtlName + "/Sim" +  str(index) + ".out"
			outputfileZe = outputPath + "/" + eqtlName + "/Sim" +  str(index) + "_eqtl.Z";
			outputfileZg = outputPath + "/" + eqtlName + "/Sim" +  str(index) + "_gwas.Z";
			saveSimData(gtexSNPIds, simZeQTL, simZGWAS, outputfileZe, outputfileZg);
			#save to file and call eCAVIAR
		
			script += ecaviarPath + "/eCAVIAR " + \
			 " -o " + outputfile + \
			 " -z " + outputfileZe + \
			 " -z " + outputfileZg + \
			 " -l " + ldFile + \
			 " -l " + ldFile + \
			 " -c 3 " + \
			 " -f 1 -r 0.95 \n";
	
		while( maxSubmitReached(400) ) :
                	time.sleep(1);	
		scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
		scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="eCAVIAR", logfile=logfile, errfile=errfile, ecaviarPath=ecaviarPath, slots=1))
		scriptFILEHandler.close();
		subprocess.call('qsub ' + scriptfile + '.qsub', shell=True)
print "DONE WE SUBMIT ALL";
