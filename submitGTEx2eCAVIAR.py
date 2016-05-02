import os
import glob
import time
import random
import sys, traceback
import subprocess
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
#$ -l h_data=3G
#$ -l h_rt=1:00:00
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
module load python
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""
gwasName = "FP";
sigFile = "../GTEx_GWAS/sign_snps_FP_MAGIC_FastingProinsulin.txt";
gtexPath = "/u/project/eeskin/zarlab/abzhu/data/45282/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2015-01-12/eqtl_data/";
gwasFile  = "~/project/data/GTEx_GWAS/Final/MAGIC/FP_MAGIC_FastingProinsulin.hg19.bed.final"
tissueNames = glob.glob(gtexPath + "GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/*.eqtl");
tissueNames = [tissueName.replace(gtexPath+"GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/", '').replace('_Analysis.cis.eqtl','') for tissueName in tissueNames];
print tissueNames;

outputPath = os.getcwd();
logFolder = "/tmp/log/"
print outputPath;
for tissueName in tissueNames:
	#make the bash file to submit
        print tissueName;
	scriptfile = outputPath + "/" + logFolder + "/command_line_" + tissueName;
	logfile    = outputPath + "/" + logFolder + "/command_line_" + tissueName + ".log";	
	errfile    = outputPath + "/" + logFolder + "/command_line_" + tissueName + ".err"; 
	if not os.path.exists(outputPath + "/tmp") :
		os.makedirs(outputPath + "/tmp");
	if not os.path.exists(outputPath + "/tmp/peak/") :
                os.makedirs(outputPath + "/tmp/peak/");
	if not os.path.exists(outputPath + "/tmp/in/") :
                os.makedirs(outputPath + "/tmp/in/");
	if not os.path.exists(outputPath + "/" + logFolder) :
		os.makedirs(outputPath + "/" + logFolder);

	script  = "python " + outputPath+ "/peakFinder.py -n " + gwasName + \
		  " -g " + gtexPath + \
		  " -f " + gwasFile + \
		  " -s " + sigFile  + \
		  " -o " + outputPath + "/tmp/" + \
		  " -t " + tissueName;
	print script;
	open(scriptfile+'.qsub', 'wb').write(TEMPLATE_SERIAL.format(script=script, name="GTEx2eCAVIAR", logfile=logfile, errfile=errfile, slots=1))
	subprocess.call('qsub -cwd -l highp,h_rt=1:00:00,h_data=3G ' + scriptfile + '.qsub', shell=True)

