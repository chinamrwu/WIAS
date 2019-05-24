import os
import sys
import time
from datetime import datetime 
startT=time.time()

step01=("OpenSwathWorkflow -min_upper_edge_dist 1 -mz_extraction_window 0.01 "
"-rt_extraction_window 600 "
"-extra_rt_extraction_window 100 "
"-min_rsq 0.8 -min_coverage 0.7 "
"-Scoring:stop_report_after_feature 5 "
"-tr /work/data/plasma20180822.TraML "
"-tr_irt /work/data/iRT.plasma.TraML "
"-threads 18 "
"-readOptions workingInMemory "
"-tempDirectory /tmp "
"-in /work/mzXML/%s "
"-force "
"-out_tsv /work/output/remedy/openswath/%s ")


step02=("mProphetScoreSelector.sh"
" /work/output/remedy/openswath/%s xx_swath_prelim_score"
" library_corr yseries_score xcorr_coelution_weighted"
" massdev_score norm_rt_score library_rmsd bseries_score"
" intensity_score xcorr_coelution log_sn_score isotope_overlap_score"
" massdev_score_weighted xcorr_shape_weighted isotope_correlation_score xcorr_shape ")

step03 = "pyprophet --ignore.invalid_score_columns --target.dir=/work/output/remedy/pyprophet --xeval.num_iter=10 --d_score.cutoff=0.01 /work/output/remedy/openswath/%s "





dsFiles=os.listdir("/work/output/pyprophet")
dsFiles=[f1 for f1 in dsFiles if f1.endswith("dscore_filtered.csv")]
dsFiles=map(lambda x:x.split("_with_descore")[0],dsFiles)
mzXMLs=[f1 for f1 in os.listdir("/work/mzXML") if f1.endswith(".mzXML")]
mzXMLs=map(lambda x:x.split(".")[0],mzXMLs)
todo=list(set(mzXMLs)-set(dsFiles))
todo=map(lambda x: x+".mzXML",todo)
#for e in todo:
#	print(e)
print("Total %d files failed running OpenSWATH " len(todo))


params=[int(sys.argv[1])-1,int(sys.argv[2])]
frm = min(params)
to  = max(params)


if frm<0:
	frm=0
if to<0:
	to=0
if to>len(todo):
	to=len(todo)
tasks=files[frm:to]
for task in tasks:                
	fname=task
	sname=fname.split(".")[0]+".tsv"
	try:
		os.system(step01 % (fname,sname))
		os.system(step02 % sname)
		os.system(step03 % sname)
		print(task)
	except:
		print("***************** Exception :"+task+"  *************************** ") 
print("========================================================== DONE ========================================================================")
endT=time.time()
duration=(endT-startT)
print("\n\n                                Time used: %d hour(s) %d minute(s) %d seconds" % (duration//3600,(duration%3600) // 60 ,(duration % 3600) % 60))
os.system("exit")

