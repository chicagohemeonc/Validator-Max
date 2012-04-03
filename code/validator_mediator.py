#!/usr/bin/env python

import sys, os, commands, string, time, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')

if sys.argv[2]=="false": #Input_data1_source=default input
    DAT_FILE = galaxyhome + "/tools/proteonics/" + sys.argv[3]    #for VM Galaxy
#    DAT_FILE = "/meida/Work/galaxy-proteonics/tools/proteonics/" + sys.argv[3]    #for test
else:
    DAT_FILE = sys.argv[3]   #Input_data1_source=uploaded input

if sys.argv[5]=="none": #Input_data2_source=none
    CPAS_FILE = ""    
elif sys.argv[5]=="default":
#    CPAS_FILE = "/meida/Work/galaxy-proteonics/tools/proteonics/" + sys.argv[6]    #for test
   CPAS_FILE = galaxyhome + "/tools/proteonics/" + sys.argv[6]    #for VM Galaxy
else:
    CPAS_FILE = sys.argv[6]   #Input_data2_source=uploaded input


PEAK_CUTOFF = sys.argv[8]
MS1_CUTOFF = sys.argv[10]
MS2_CUTOFF = sys.argv[12]
MAX_PEAK_DIFF = sys.argv[14]

htmlout = sys.argv[16]  #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66.dat

OUTPUT_FOLDER = htmlout.split('.')[0]+'_files'   #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66_files

if not os.path.exists(OUTPUT_FOLDER):
	os.makedirs(OUTPUT_FOLDER)
else:
	pass

print "python","validator_cli.py","--input-file1",DAT_FILE,"--input-file2",CPAS_FILE,"--peak-cutoff",PEAK_CUTOFF,"--ms1-cutoff", MS1_CUTOFF,"--ms2-cutoff", MS2_CUTOFF,"--max-peak-diff", MAX_PEAK_DIFF,"--output-folder", OUTPUT_FOLDER


'''
usage: validator_cli.py [-h] [--input-file1 DAT_FILE]
	[--input-file2 CPAS_FILE]
	[--peak-cutoff PEAK_CUTOFF]
	[--ms1-cutoff MS1_CUTOFF]
	[--ms2-cutoff MS2_CUTOFF]
	[--max-peak-diff MAX_PEAK_DIFF]
	[--output-folder OUTPUT_FOLDER]
'''


galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlattr = """Galaxy tool %s run at %s</b><br/>"""
galhtmlpostfix = """</div></body></html>\n"""


def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def htmloutput(htmlout,outputfolder):
	rstyle="""<style type="text/css">
        tr.d0 td {background-color: oldlace; color: black;}
        tr.d1 td {background-color: aliceblue; color: black;}
        </style>"""    
        res = [rstyle,]
        res.append(galhtmlprefix % os.path.basename(sys.argv[0]))   
        res.append(galhtmlattr % ('validator',timenow()))
        flist = [x for x in os.listdir(outputfolder) if not x.startswith('.')] 
        if len(flist) > 0:
            res.append('<b>The following output files were created (click the filename to view/download a copy):</b><hr/>')
            res.append('<table>\n')
            for i,f in enumerate(flist):
                fn = os.path.split(f)[-1]
                res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (fn,fn))
            res.append('</table><p/>\n') 

        res.append(galhtmlpostfix) 
        outf = open(htmlout,'w')
        outf.write(''.join(res))   
        outf.write('\n')
        outf.close()


os.chdir(galaxyhome + "/tools/proteonics/")

if CPAS_FILE == '':
	subprocess.call(["python","validator_cli.py","--input-file1",DAT_FILE,"--peak-cutoff",PEAK_CUTOFF,"--ms1-cutoff", MS1_CUTOFF,"--ms2-cutoff", MS2_CUTOFF,"--max-peak-diff", MAX_PEAK_DIFF,"--output-folder", OUTPUT_FOLDER])
else:
	subprocess.call(["python","validator_cli.py","--input-file1",DAT_FILE,"--input-file2", CPAS_FILE, "--peak-cutoff",PEAK_CUTOFF,"--ms1-cutoff", MS1_CUTOFF,"--ms2-cutoff", MS2_CUTOFF,"--max-peak-diff", MAX_PEAK_DIFF,"--output-folder", OUTPUT_FOLDER])



htmloutput(htmlout,OUTPUT_FOLDER)


