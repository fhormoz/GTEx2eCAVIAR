import os
import time
import random
import glob
import numpy as npy
import sys, traceback
import subprocess
from myconfig import *
from subprocess import Popen, PIPE


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
#$ -l h_data=2G
#$ -l h_rt=4:00:00
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
outputPath  = currentPath + "/output/" + gwasName + "/";

if not os.path.exists(currentPath):      
        os.makedirs(currentPath);
if not os.path.exists(currentPath+"/output/"):
        os.makedirs(currentPath+"/output/");
if not os.path.exists(outputPath):
	os.makedirs(outputPath);
if not os.path.exists(outputPath+"/log/"):      
        os.makedirs(outputPath+"/log/");

for tissue in glob.glob(currentPath+"/tmp/in/" + gwasName + "/*") :
	while( maxSubmitReached(400) ) :
		time.sleep(4);
	#make the bash file to submit
	for eqtlFile in glob.glob(tissue + "/*.eqtl"):
		tissueName = tissue.replace(currentPath+"/tmp/in/" + gwasName + "/", '');
		while( maxSubmitReached(400) ) :
                	time.sleep(4);
		eqtlName = eqtlFile.replace(currentPath + "/tmp/in/" + gwasName + "/" + tissueName + "/",'').replace('.eqtl','');	
		scriptfile = outputPath + "/log/command_" + eqtlName;
		logfile    = outputPath + "/log/command_" + eqtlName + ".log";
		errfile    = outputPath + "/log/command_" + eqtlName + ".err";
		
		outputfile = outputPath + "/" + eqtlName + ".out";
		ldfile     = eqtlFile.replace('.eqtl', '.LD');
		zfile      = eqtlFile.replace('.eqtl', '.gwas'); 
		eQTLdata   = npy.genfromtxt(eqtlFile, usecols=(1,));
		zGWASdata  = npy.genfromtxt(zfile, usecols=(1,));
		if(npy.ndim(eQTLdata) ==0):
			continue;
		if( max(abs(eQTLdata)) < 4 ):
			continue;
		if( max(abs(zGWASdata)) < 4 ):
                        continue;
		if os.path.exists(outputfile + "_col"):
			continue;
		print outputfile;
		
		script = ecaviarPath + "/eCAVIAR " + \
                 " -o " + outputfile + \
                 " -z " + eqtlFile + \
                 " -z " + zfile + \
                 " -l " + ldfile + \
		 " -l " + ldfile + \
                 " -c 3 " + \
                 " -f 1 -r 0.95";
		print script;
		scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
		scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="eCAVIAR", logfile=logfile, errfile=errfile, ecaviarPath=ecaviarPath, slots=1))
		scriptFILEHandler.close();
		subprocess.call('qsub -v PYTHONPATH,LIBRARY_PATH ' + scriptfile + '.qsub', shell=True)
	
