import sys
import glob
import numpy as np
import scipy.stats
import pdb
import warnings

posToRS = {}
for l in open("/project/voight_MR/ksiewert/LDSC_LipidMVP/LDSC_LipidMVP_imputation0.8/MVP_HDL_plinkform_withrs.txt","r"):
	sl = l.split(" ")
	posToRS["chr"+sl[0]+":"+sl[1]] = sl[2]


output = open("/project/voight_MR/ksiewert/MenarcheClumps/formMenarcheSumstats.txt",'w')
output.write(" ".join(["CHR","BP","SNP","A1","A2","BETA","SE","P","\n"]))
for l in open("/project/voight_datasets/GWAS/51_menarche/Menarche_1KG_NatGen2017_WebsiteUpload.txt",'r'):
	sl = l.split("\t")
	if sl[0] in posToRS:
		RS = posToRS[sl[0]]
		p = float(sl[4])
		B = float(sl[3])
		z = scipy.stats.norm.ppf(p/2.)
		A1 = sl[1]
		A2 = sl[2]
		with warnings.catch_warnings():
			warnings.filterwarnings('error')
			try:
				se = str(abs(B)/abs(z))
			except Warning:
				se = "NaN" 
		ch = sl[0].split(":")[0].replace("chr","")
		bp = sl[0].split(":")[1]
		output.write(" ".join([ch,bp,RS,A1,A2,str(B),se,str(p)])+"\n")

output.close()
