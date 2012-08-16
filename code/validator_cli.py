####
#Copyright 2011 Samuel Volchenboum, Jonathan Goya, Gene Selkov, Chaim Kirby, 
#
#This file is part of Validator.
#
#Validator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Validator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with Validator.  If not, see <http://www.gnu.org/licenses/>.
#####
#!/usr/bin/env python

import sys
import cPickle
import time
from PeptideFragmentSingleton import PeptideFragment 
import bisect
import os
import itertools # built in to 2.6
from xlwt import Workbook,Style,easyxf # get these from http://www.python-excel.org/
from math import trunc
import argparse
import subprocess
try:
	from padnums import pprint_table
except:
	pass

path =  os.path.abspath(os.path.dirname(__file__))+"/../"

parser = argparse.ArgumentParser()

parser.add_argument('--input-file1', 
		    action='store', 
                    dest='dat_file', 
                    default=path+'data/sample_data/F004642_sec_12.dat', 
                    help='A string that points to the location of input file')

parser.add_argument('--input-file2',
                    action='store',
                    dest='cpas_file',
                    default='',
                    help='A string that points to the location of cpas input file [Required for CROSSREF and QUANT to run]')


parser.add_argument('--peak-cutoff',
                    action='store',
                    dest='peak_cutoff',
                    default=0.01, 
                    help='Peak Cutoff value. Default is 0.01')

parser.add_argument('--ms1-cutoff',
                    action='store',    
                    dest='ms1_cutoff',
                    default=10,
                    help='MS1 Cutoff value. Default is 10')

parser.add_argument('--ms2-cutoff',
                    action='store',
                    dest='ms2_cutoff',
                    default=1000,
                    help='MS2 Cutoff value. Default is 1000')

parser.add_argument('--max-peak-diff',
                    action='store',
                    dest='max_peak_diff',
                    default=1000,
                    help='Max Peak Difference value. Default is 0.25')

parser.add_argument('--output-folder',
                    action='store',
                    dest='output_folder',
                    default=path+'output_files/sample_data/F004642_sec_12/',
                    help='Output folder location')


results = parser.parse_args()

global screen_print
screen_print2 = False
screen_print = False


#global data_folder
#data_folder = sys.argv[1]
#data_folder = 'sample_data'

global input_file
input_file = results.dat_file

global cpas_file
cpas_file = results.cpas_file

global output_folder
output_folder = results.output_folder

global csv_name
global txt_name
global xls_name
csv_name = "Val_1_pairs.csv"
txt_name = "information_file.txt"
xls_name = "validator_3.xls"

global output_csv
global output_txt
global output_xls
output_csv = output_folder+"/"+csv_name
output_txt = output_folder+"/"+txt_name
output_xls = output_folder+"/"+xls_name

global peak_cutoff
#
peak_cutoff = results.peak_cutoff

global MS1_cutoff
#
MS1_cutoff = results.ms1_cutoff

global MS2_cutoff
#
MS2_cutoff = results.ms2_cutoff
	
global max_peak_difference
#
max_peak_difference = results.max_peak_diff

global proton_mass
proton_mass = 1.00727646688
neutron_mass = 1.0086649156
first_isotopologue = 1.0033548378 # http://en.wikipedia.org/wiki/Isotopes_of_carbon
second_isotopologue  = 2.003241989
ppm_cutoff = int(MS1_cutoff) # ppm
scan_width = 100
ms2_ppm_tolerance = int(MS2_cutoff)
min_ion_percent_intensity = float(peak_cutoff)# mininum to keep
max_ion_difference = float(max_peak_difference) # max allowable difference between ion percentages to keep

##########################
#
#  Stanard functions
#
##########################

def exit():
	sys.exit()

def load_pickle_data(filename):
	if screen_print:
		print "Loading saved data...",filename
	newfile = open(filename)
	data = cPickle.load(newfile)
	newfile.close()
	return data

def print_timing(func):
	def wrapper(*arg):
		t1 = time.time()
		res = func(*arg)
		t2 = time.time()
		print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
		return res
	return wrapper

def writefile(output,filename):
	newfile = open(filename, "w")
	for a in output:
		newfile.write(a+"\n")
	newfile.close()
	return

def pickle_data(data, filename):
	newfile = open(filename,"w")
	cPickle.dump(data,newfile)
	newfile.close()
	return
	
def readfile(name_of_file,cr = 1):
	if screen_print:
		print "Loading "+name_of_file+"..."
	file_data = open(name_of_file)
	file_contents = file_data.readlines()
	output = []
	for i,line in enumerate(file_contents):
			if (i <> len(file_contents)-1): # no CR on last line, apparently
					output.append(line[:(-1* cr)]) # strip carriage returns
			else:
					output.append(line)
	if screen_print:
		print "...done"
	return output
	
def make_data_dictionary(data): # creates dictionary where keys are section names
	if screen_print:
		print "Creating data dictionary..."
	content = {}
	flag = 0
	for line in data:
			if "Content-Type: application" in line: #new section
					flag = 1
					section_name = line.split('"')[1] #key - eg. peptides, parameters
					if "query" in section_name:
							flag = section_name[5:] # query number
							section_name = "query"+flag
					continue
			elif flag > 0:
					if flag == 1: # new section, not a query
						content.setdefault(section_name,[]).append(line)
					else:						
						if line == "" or len(line) == 1:
							content.setdefault(section_name,[]).append("query"+str(flag))
						else:
							content.setdefault(section_name,[]).append(line)
	return content

def make_query_dictionary(mass_data,peptide_data,protein_data):

	def make_peptide_dictionary(data):
		if screen_print:
			print "Creating peptide dictionary..."
		output = {}
		for line in data[1:-1]:

			if "terms" not in line and "=-1" not in line and "primary" not in line and "subst" not in line:
				key = int(line.split('_')[0][1:]) # query number
				pep = line.split(';')[0].split(',')
				protein_info = line.split(';')[1]
				protein_info = protein_info.replace(',','%')
				pep.append(protein_info)
				#print pep
				val = peptide_line_object(pep)
				#print val.list
				output.setdefault(key,[]).append(val)
		return output

	def make_MS2(data):
		# keys are query numbers
		# many will have the same scan numbers
		if screen_print:
			print "Creating protein dictionary..."
			
		output = {}
		for i in data:
			if len(i) == 11: # this is a consistency check to make sure that all the MS2 scan data is there. If not, don't bother.				
				key = int(i[0][5:])
				good_queries.append(key)	
				val = i[1:10]
				output.setdefault(key,[]).append(val)
		return output
	good_queries = [] 
	output = {}
	peptide_dictionary = make_peptide_dictionary(peptide_data)
	peptide_queries = peptide_dictionary.keys()
	peptide_queries.sort()
	protein_dictionary = make_MS2(protein_data)
	protein_queries = protein_dictionary.keys()
	for i,line in enumerate(mass_data):
		if 'qexp' in line:
			query_number = int(line.split('=')[0][4:])
			if query_number in good_queries:
				mass = float(line.split('=')[1].split(',')[0])
				charge = float(line.split(',')[1][0])
				qmass = mass_data[i-1]
				qmatch = mass_data[i+1]
				qplughole = mass_data[i+2]
				precursor_mass = mass * charge - 1.007825*charge
				key = str(query_number)
				val = [[query_number,mass,charge,precursor_mass,qmass,line,qmatch,qplughole]]
				try:
					val.append(peptide_dictionary[query_number])
				except KeyError:
					val.append([])
				try:
					val.append(protein_dictionary[query_number][0])
				except KeyError:
					val.append(["No proteins"])
				#print key
				val = q_dict_object(val)
				output.setdefault(key,[]).append(val)
	return output

def ppm(m1,m2):
	"""calculated parts per million (ppm) difference betwen two masses"""
	return float(abs(m1-m2))/float(min(m1,m2)) * 1000000.0	

def calc_ppm(mass1, mass2):
	mass1 = float(mass1)
	mass2 = float(mass2)
	ppm = 1000000 * (mass2-mass1)/mass1
	#print "PPM",ppm
	return ppm

def create_peptide_dictionary(dd):
	keys = dd.keys()
	new_peptide_dictionary = {}
	for key in keys:
		element = dd[key][0]
		key_ = element[0][0]
		new_peptides  = element[1]
		for new_peptide in new_peptides:
			new_peptide_dictionary.setdefault(key_,[]).append(new_peptide)
	return new_peptide_dictionary

def unique(s):
    n = len(s)
    if n == 0:
        return []
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

def same_ion(ion1,ion2):
	ppm = calc_ppm(ion1,ion2)
	if abs(ppm) < ms2_ppm_tolerance:
		return True
	else:
		return False

def scan_number(dictionary_element):
	'''Returns scan number from a dictionary element'''
	#return int(dictionary_element[2][0].split('FinneganScanNumber')[1].split('%20')[1]) # might have to use this line depending on Mascot version
	return int(dictionary_element[2][0].split('%2e')[1])

def extract(d,keys):
	
	# used to create a smaller dictionary subset for testing
	return dict((k, d[k]) for k in keys if k in d)

def average(somelist):
	if len(somelist)>0:
		return float(sum(somelist))/float(len(somelist))
	else:
		return 0

def flatten(x):
    """flatten(sequence) -> list"""
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def av_ppm(somelist, i):
	if i == "B":
		subt = 0
	else:
		subt = isotope_masses[i]
	if len(somelist) > 0:
		return sum(ppm(a[1][0]-subt,a[0][0]) for a in somelist) / float(len(somelist))
	else:
		return 0

class q_dict_object(list):
	def __init__(self, data):
		list.__init__(self, data)
		self.charge = int(data[0][2])
		#self.precursor_MW = float(data[0][3])
		self.precursor_mz = float(data[0][1])
		self.precursor_MW = (self.precursor_mz * self.charge) - (self.charge * proton_mass)
		self.peptide_lines = data[1]
		self.peptides = [a.peptide for a in data[1]]
		self.proteins = [a.protein for a in data[1]]
		self.scores = [a.score for a in data[1]]
		if len(data[1]) > 0:
			self.top_peptide_line = data[1][0]
			self.top_peptide = data[1][0].peptide
			self.top_PME = data[1][0].PME
			self.top_mods = data[1][0].mods
			self.top_mod = data[1][0].mod
			self.top_score = data[1][0].score
			self.top_protein = data[1][0].protein
			self.top_peaks_used = data[1][0].peaks_used
			self.top_calc_mz = data[1][0].calc_mz
		else:
			self.top_peptide = "None"
			self.top_score = 0
			self.top_protein = "None"
			self.top_mods = '0'
			self.top_mod = '0'
			self.top_PME = 0
			self.top_calc_mz = 0
			
		self.MS2_data = data[2]
		self.MS2_ions = [[float(b) for b in a.split(':')] for a in data[2][8][6:].split(',')] #.sort(lambda x,y:cmp(x[0],y[0]))
		self.scan_number = scan_number(data)
class peptide_line_object(list):
	def __init__(self,data):
		list.__init__(self, data)
		self.calc_mz = float(data[1])
		self.mods = data[6]
		self.mod = int(data[6][-2])
		self.score = float(data[7])
		self.peptide = data[4]
		self.PME = float(data[2])
		self.score = float(data[7])
		#print self.peptide
		#print data[11]
		self.protein = [a.split(':')[0].replace('"','') for a in data[11].split('%')]
		#self.protein = [a.split('"')[1].split('|')[1] for a in data[11].split('%')]
		#self.protein = data[11].split('"')[1].split('|')[1] # yeast
		self.peaks_used = int(data[5])


def print_above_and_below_35(queries,Q):
	q_above = [q for q in queries if len(Q[q][0].peptides)>0 and Q[q][0].top_score >= 35]
	q_below = [q for q in queries if len(Q[q][0].peptides)>0 and Q[q][0].top_score < 35]
	print "There are ",len(q_above),"queries above 35 and",len(q_below),"queries below 35"
	pep_above = unique([Q[q][0].top_peptide for q in q_above])
	pep_below = unique([Q[q][0].top_peptide for q in q_below])
	print "There are ",len(pep_above),"peptides above 35 and",len(pep_below),"peptides below 35"
	return


def mods_same(mods1,mods2):
	'''Takes two sets of mods and determines if is possible that one is a light and other is a heavy isotope'''
	if mods1 == mods2:
		return False
	elif len(mods1) == len(mods2) and (( len([a for a in mods_to_check if str(a) in mods1])>0 and len([a for a in mods_to_check if str(a) in mods2]) == 0 ) or  ( len([a for a in mods_to_check if str(a) in mods2])>0 and len([a for a in mods_to_check if str(a) in mods1]) == 0 )) :
		for i,m1 in enumerate(mods1):
			m1 = int(m1)
			m2 = int(mods2[i])
			# print m1, m2, int(m1) in mods_to_check
			if (m1<>m2 and not((m1 == 0 and m2 in mods_to_check) or (m1 in mods_to_check and m2 == 0))) or (m1 in mods_to_check and m2 in mods_to_check):
				return False
		return True	
	return False

def mods_consistent(pep,mods,kind):

	M = ['K','R']
	
	for l in zip(pep,mods):
		p = l[0]
		m = int(l[1])
		if p in M and kind == 'L' and m <> 0:
			return False
		elif p in M and kind == 'H' and m <> mod_dict[p]:
			return False
		elif p not in M and m in mods_to_check:
			return False
	return True
		
def validator_1(Q):
	if screen_print:
		print
		print "##################################"
		print "#########   VALIDATOR 1 ##########"
		print "##################################"
		print
	queries = Q.keys()
	queries.sort(lambda x,y: cmp(Q[x][0].scan_number,Q[y][0].scan_number))
	scan_numbers = [Q[a][0].scan_number for a in queries]
	D = zip(scan_numbers,queries)
	pairs = []

	for i,line in enumerate(D):
		q1 = line[1]
		c1 = Q[q1][0].charge
		pep1 = Q[q1][0].top_peptide
		s1 = line[0]
		###
		### Find the scan window 
		###
		low = max(s1-scan_width,0)
		while low not in scan_numbers and low <=s1:
			low += 1
		high = min(s1+scan_width, scan_numbers[-1])
		while high not in scan_numbers and high >= s1:
			high -= 1
		i1 = scan_numbers.index(low)
		i2 = scan_numbers.index(high)
		R2 = D[i1:i2]
		#
		#
		#print mods_to_check		
		for line2 in R2: # iterate through the scans
			q2 = line2[1]
			mz1 = Q[q1][0].precursor_mz
			mz2 = Q[q2][0].precursor_mz
			pep2 = Q[q2][0].top_peptide
			c2 = Q[q2][0].charge
			s2 = line2[0]
	
			if q1 <> q2 and pep1 == pep2 and c1 == c2 and mz1<>mz2:
				if mz1<mz2:
					pepL = pep1
					pepH = pep2
					modsL = Q[q1][0].top_mods[1:-1]
					modsH = Q[q2][0].top_mods[1:-1]
					qL = q1
					qH = q2
				else:
					pepL = pep2
					pepH = pep1
					modsL = Q[q2][0].top_mods[1:-1]
					modsH = Q[q1][0].top_mods[1:-1]
					qL = q2
					qH = q1
				
				if modsL <> modsH and len(pepL) == len(modsL) and len(pepH) == len(modsH) and mods_consistent(pepL,modsL,'L') and mods_consistent(pepH,modsH,'H') and [qL,qH] not in pairs: 
					pairs.append([qL,qH])

	pairs.sort(lambda a,b:cmp(int(Q[a[0]][0].scan_number),int(Q[b[0]][0].scan_number)))
	output = []
	line1 = "Query 1, Query 2, Peptide, Mods 1, Mods 2, Scan 1, Scan 2, Score 1, Score 2, PME, ? Isotopologue, Proteins"
	output.append(line1)
	for line in pairs:		
		q1 = line[0]
		q2 = line[1]		
		mods1 = Q[q1][0].top_mods[1:-1]
		mods2 = Q[q2][0].top_mods[1:-1]
		dif_modsL = [a[0] for a in zip(mods1,mods2) if a[0] <> a[1]] # find the differences between the two 
		dif_modsH = [a[1] for a in zip(mods1,mods2) if a[0] <> a[1]] # find the differences between the two 
		subL = sum([modifications[int(mod)][1] for mod in dif_modsL])
		subH = sum([modifications[int(mod)][1] for mod in dif_modsH])
		sub = subH-subL
		ME = (Q[q2][0].precursor_MW - sub) - Q[q1][0].precursor_MW 
		PME = ME - neutron_mass * (ME > 0.5) + neutron_mass * ( ME < -0.5 )
		if ME > 0.5 or ME < -0.5:
			add = "**"
		else:
			add = ""		
		ppm_calc = PME/Q[q1][0].precursor_MW *1000000
		proteins = Q[q1][0].top_protein
		prots = ",".join(proteins)		
		score1 = Q[q1][0].top_score
		score2 = Q[q2][0].top_score		
		newline =  [q1,q2,Q[q1][0].top_peptide,"'"+str(mods1),"'"+str(mods2),Q[q1][0].scan_number,Q[q2][0].scan_number,score1, score2,abs(ppm_calc),add, prots]
		newline = ",".join([str(a) for a in newline])
		output.append(newline)
	writefile(output,output_csv)	
	return

def tally_match_score(i):
	'''This takes light and heavy ions and calculates how many are found in both'''
	light = i[0]
	heavy = i[1]
	ranges = [[0,4],[5,9]]
	tallies = [0,0]
	for x in zip(light, heavy):
		for z in [0,1]:
			for y in range(ranges[z][0],ranges[z][1]):
				if type(x[0][y]) == float and type(x[1][y]) == float:
					tallies[z] += 1			
	return tallies


def writeexcel(output,filename):
	# this section written to write Validator 3 data out to an excel workbook
	# all documentation from http://www.python-excel.org/ and examples therein

	header = [a[0] for a in output[0]]
	
	widths = [a[1] for a in output[0]]
	# these are the row titles for the fragmentaion page
	# they have to match up with the column headings
	header_1 = [a[3] for a in output[0] if a[2] == 0]

	# get list of peptides that are in good matches
	# this allows us to keep track of which queries are part of good matches
	peptide_L_index = header.index("Peptide L")
	peptide_H_index = header.index("Peptide H")
	isotope_match_L_index = header.index("Isotope match L")
	isotope_match_H_index = header.index("Isotope match H")
	isotope_and_mods_match_L_index = header.index("Isotope + Mods match L")
	isotope_and_mods_match_H_index = header.index("Isotope + Mods match H")
	mods_match_L_index = header.index("Mods match L")
	mods_match_H_index = header.index("Mods match H")
	match_index = header.index("Match?")
	pep_found_L_index = header.index("Pep found L")
	pep_found_H_index = header.index("Pep found H")
	peptides_in_match = []
	for line in output[1:]:
		d = line[0]
		if d[match_index] == "True":
			peptides_in_match += [d[peptide_L_index], d[peptide_H_index]]
	peptides_in_match = unique(peptides_in_match)

	'''colors
	aqua black blue blue_gray bright_green brown
	coral cyan_ega dark_blue dark_blue_ega dark_green dark_green_ega dark_purple dark_red
	dark_red_ega dark_teal dark_yellow gold
	gray_ega gray25 gray40 gray50 gray80 green ice_blue indigo ivory lavender
	light_blue light_green light_orange light_turquoise light_yellow lime magenta_ega ocean_blue olive_ega olive_green orange pale_blue periwinkle pink
	plum purple_ega red rose sea_green silver_ega sky_blue tan teal teal_ega turquoise violet white yellow
	'''
	wb = Workbook()
	ws_pair_data = wb.add_sheet('pair_data',cell_overwrite_ok=True)
	ws_pair_data.panes_frozen = True
	ws_pair_data.remove_splits = True
	ws_pair_data.horz_split_pos = 1
	ws_pair_data.horz_split_first_visible = 2
	fragmentation_data_sheets = 1
	ws_fragmentation_data = wb.add_sheet('fragmentation_data_'+str(fragmentation_data_sheets),cell_overwrite_ok=True)
	# excel styles
	style_none = easyxf('font: name Courier;')
	style = easyxf('font: name Courier; align: horizontal center;')
	dec_style = easyxf('font: name Courier;', num_format_str='0.000')
	dec_style_highlight = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour yellow;', num_format_str='0.000')
	dec_style_no_highlight = easyxf('font: name Courier; align: horizontal center;', num_format_str='0.000')
	color_row_blue = easyxf('font: name Courier; pattern: pattern solid, fore_colour blue')
	style_light_green = easyxf('font: name Courier; pattern: pattern solid, fore_colour light_green')
	style_center_light_green = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour light_green')
	style_center_light_red = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour red')
	# sheet 2 formatting
	color_y = "light_yellow"
	style_colored_y = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour ' + color_y + ';')
	decimal_style_colored_y = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour ' + color_y + ';', num_format_str='0.000')
	color_b = "light_turquoise"
	style_colored_b = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour ' + color_b + ';')
	decimal_style_colored_b = easyxf('font: name Courier; align: horizontal center; pattern: pattern solid, fore_colour ' + color_b + ';', num_format_str='0.000')	
	decimal_style_no_color = easyxf('font: name Courier; align: horizontal center;', num_format_str='0.000')
	style_no_color = easyxf('font: name Courier; align: horizontal center;')
	# border formatting
	border_tl = easyxf('border: left thick, top thick')
	border_t = easyxf('border: top thick')
	border_tr = easyxf('border: right thick, top thick')
	border_r = easyxf('border: right thick')
	border_br = easyxf('border: right thick, bottom thick')
	border_b = easyxf('border: bottom thick')
	border_bl = easyxf('border: left thick, bottom thick')
	border_l = easyxf('border: left thick')

	row1 = 0	
	for col,item in enumerate(header):
		ws_pair_data.row(row1).write(col,item,style)
		ws_pair_data.col(col).width = 256*widths[col]
	row1 += 1
	row2 = 0
	ws_fragmentation_data.col(0).width = 256*20	
	for col in [0,1,2,4,15]:
		ws_fragmentation_data.col(col).width = 256*30
	ws_fragmentation_data.col(3).width = 256*5

	#########################
	#                       #
	#  MAIN OUTPUT SECTION  #
	#                       #
	#########################

	for line in output[1:]:
		data = line[0]
		for i,l in enumerate(data):
			if l == "False":
				data[i] = ""
			elif l == "True":
				data[i] = "T"
		ion_set = line[1] # will have lengh of 0, 1, 2
		# if the peptides are the same, ion_set will have length 1
		# if the peptides are different but the mods match, then ion_set will have length of 1 or 2, depending on if one or both match			
		ions = []
		consecutive = []
		scores = []
		tallies = []
		pc = []
		c = ""

		for i in ion_set:
			I = [a[0] for a in i]	
			C = [a[1] for a in i] 
			S = [a[2] for a in i]
			P = [a[3] for a in i]
			T = tally_match_score(I)
			ions.append(I)
			consecutive.append(C)
			scores.append(S)
			tallies.append(T)
			pc.append(P) # percent coverage
			
		if data[match_index] == "T":
			data.append(scores[0][0][0])
			data.append(scores[0][0][1])
			data.append(sum(scores[0][0]))
			data.append(consecutive[0][0])
			data.append(pc[0][0][0])
			data.append(pc[0][0][1])
			data.append(sum(pc[0][0]))
			data.append(scores[0][1][0])
			data.append(scores[0][1][1])
			data.append(sum(scores[0][1]))
			data.append(consecutive[0][1])
			data.append(pc[0][1][0])
			data.append(pc[0][1][1])
			data.append(sum(pc[0][1]))
			data.append(tallies[0][0])
			data.append(tallies[0][1])
			data.append(tallies[0][0] + 2 * tallies[0][1])
			data += [c] * 20
			
		elif data[isotope_and_mods_match_L_index] == "T" and data[isotope_and_mods_match_H_index] <> "T":
			data.append(scores[0][0][0])
			data.append(scores[0][0][1])
			data.append(sum(scores[0][0]))
			data.append(consecutive[0][0])
			data.append(pc[0][0][0])
			data.append(pc[0][0][1])
			data.append(sum(pc[0][0]))
			data += [c] * 10
			data.append(scores[0][1][0])
			data.append(scores[0][1][1])
			data.append(sum(scores[0][1]))
			data.append(consecutive[0][1])
			data.append(pc[0][1][0])
			data.append(pc[0][1][1])
			data.append(sum(pc[0][1]))
			data.append(tallies[0][0])
			data.append(tallies[0][1])
			data.append(tallies[0][0] + 2 * tallies[0][1])
			data += [c] * 10
			
			
		elif data[isotope_and_mods_match_L_index] <> "T" and data[isotope_and_mods_match_H_index] == "T" :
			# print "HERE 3"			
			data += [c] * 7
			data.append(scores[0][1][0])
			data.append(scores[0][1][1])
			data.append(sum(scores[0][1]))
			data.append(consecutive[0][1])
			data.append(pc[0][0][0])
			data.append(pc[0][0][1])
			data.append(sum(pc[0][0]))
			data += [c] * 13
			data.append(scores[0][0][0])
			data.append(scores[0][0][1])
			data.append(sum(scores[0][0]))
			data.append(consecutive[0][0])
			data.append(pc[0][1][0])
			data.append(pc[0][1][1])
			data.append(sum(pc[0][1]))						
			data.append(tallies[0][0])
			data.append(tallies[0][1])
			data.append(tallies[0][0] + 2 * tallies[0][1])

		elif data[isotope_and_mods_match_L_index] == "T" and data[isotope_and_mods_match_H_index] == "T":
			data.append(scores[0][0][0])
			data.append(scores[0][0][1])
			data.append(sum(scores[0][0]))
			data.append(consecutive[0][0])
			data.append(pc[0][0][0])
			data.append(pc[0][0][1])
			data.append(sum(pc[0][0]))
			data.append(scores[1][1][0])
			data.append(scores[1][1][1])
			data.append(sum(scores[1][1]))
			data.append(consecutive[1][1])
			data.append(pc[0][1][0])
			data.append(pc[0][1][1])
			data.append(sum(pc[0][1]))						
			data += [c] * 3
			data.append(scores[0][1][0])
			data.append(scores[0][1][1])
			data.append(sum(scores[0][1]))
			data.append(consecutive[0][1])
			data.append(pc[1][0][0])
			data.append(pc[1][0][1])
			data.append(sum(pc[1][0]))
			data.append(tallies[0][0])
			data.append(tallies[0][1])
			data.append(tallies[0][0] + 2 * tallies[0][1])
			data.append(scores[1][0][0])
			data.append(scores[1][0][1])
			data.append(sum(scores[1][0]))
			data.append(consecutive[1][0])
			data.append(pc[1][1][0])
			data.append(pc[1][1][1])
			data.append(sum(pc[1][1]))						
			data.append(tallies[1][0])
			data.append(tallies[1][1])
			data.append(tallies[1][0] + 2 * tallies[1][1])
		else:
			data += [c] * 37
			

		#######################################################	
		#                                                     # 
		#          WRITE THE ROW IN THE DATA TABLE            #
		#                                                     # 
		#######################################################	
		if data[match_index] == "T":
			S1 = decimal_style_colored_y
			S2 = style_colored_y
		elif (data[isotope_and_mods_match_L_index] == "T" and data[pep_found_L_index] <> "T" ) or (data[isotope_and_mods_match_H_index] == "T" and data[pep_found_H_index] <> "T"):
		#elif (data[isotope_match_L_index] == "T" and data[mods_match_L_index] == "T" and data[pep_found_L_index] <> "T" ) or (data[isotope_match_H_index] == "T" and data[mods_match_H_index] == "T" and data[pep_found_H_index] <> "T"):
			S1 = decimal_style_colored_b
			S2 = style_colored_b
		else:
			S1 = decimal_style_no_color
			S2 = style_no_color
		
			
		for col,item in enumerate(header):
			if item == "Score" or item == "ABS PME":
			 	ws_pair_data.row(row1).write(col,float(data[col]),S1)
			else:
				ws_pair_data.row(row1).write(col,data[col],S2)						
		row1 += 1

		#######################################################	
		#                                                     # 
		#   WRITE THE DATA COLUMN IN THE FRAGMENTATION TABLE  #
		#                                                     # 
		#######################################################	
		if len(ions)>0:
			start = row2
			for x,item in enumerate(output[0]): # these are the header items
				if item[2] == 0: # indicates a new row
					# new row
					row2 += 1
					row_name = item[3]
					# set up styles based on whether it is a peptide match or whether the peptide is not a match but has been found before
					if row_name == "Match?" and data[x] == "T": # highlight the match
						S = [style_light_green,style_center_light_green]
					elif row_name == "Peptide" and data[x] <> data[x+1] and data[x] in peptides_in_match:
						S = [style_none,style_center_light_red]
					else:
						S = [style_none,style]
					
					ws_fragmentation_data.row(row2).write(0,row_name,S[0]) # row title
					ws_fragmentation_data.row(row2).write(1,data[x],S[1]) # the actual data

				elif row_name == "Peptide" and data[x] <> data[x-1] and data[x] in peptides_in_match:
					ws_fragmentation_data.row(row2).write(2,data[x], style_center_light_red)
				else:
					ws_fragmentation_data.row(row2).write(2,data[x], style)
		
			#######################################################	
			#                                                     # 
			#          WRITE THE FRAGMENTATION PATTERNS           #
			#                                                     # 
			#######################################################	
		
			''' Going to write either one fragmentaion table (peptides match, peptides don't match but the isotope does for one) or two (peptides don't match but both isotopes do)'''

			def write_fragmentation(r,c,i):
				length = len([a[4] for a in i[0]])
				offset =  [r, r+length+2]
				light = i[0]
				heavy = i[1]

				peptides = ["".join([a[4] for a in light]), "".join([a[4] for a in heavy])]

				for xx in range(2):
					for col,L in enumerate([peptides[xx],"B","B+2H","B-NH3","B-H2O","AA","Y","Y+2H","Y-NH3","Y-H2O"]):
						ws_fragmentation_data.row(offset[xx]).write(col+c,L,dec_style_no_highlight)
				for x,l in enumerate(light):
					for y,item1 in enumerate(l): # y = 0,1,2,3 = b_ion; y = 5,6,7,8 = y_ion
						item2 = heavy[x][y]
						if type(item1) == float and type(item2) == float: # this is a match
							ws_fragmentation_data.row(offset[0]+x+1).write(y+c+1,item1,dec_style_highlight)
							ws_fragmentation_data.row(offset[1]+x+1).write(y+c+1,item2,dec_style_highlight)
						else:
							ws_fragmentation_data.row(offset[0]+x+1).write(y+c+1,item1,dec_style_no_highlight)
							ws_fragmentation_data.row(offset[1]+x+1).write(y+c+1,item2,dec_style_no_highlight)
				return
		
			row2 = start + 1
			for x,I in enumerate(ions):	
				write_fragmentation(row2,4+11*x,I)
		
			if len(ions) > 0:
				longest = max([len(a[0]) for a in ions])
			else:
				longest = 0
			if 2 * longest < len(header_1):
				row2 += len(header_1) + 2
			else:
				row2 += 2*longest + 3
		
			for i in range(14):
					ws_fragmentation_data.row(row2).write(i,None,border_b)
			row2 += 2
			if row2 > 65000: # we've reached the limit of an .xls file, need to split into multiple sheets
				fragmentation_data_sheets += 1
				ws_fragmentation_data = wb.add_sheet('fragmentation_data_'+str(fragmentation_data_sheets),cell_overwrite_ok=True)
				ws_fragmentation_data.col(0).width = 256*20	
				for col in [0,1,2,4,15]:
					ws_fragmentation_data.col(col).width = 256*30
				ws_fragmentation_data.col(3).width = 256*5
				row2 = 0

	wb.save(filename)

	return

#@print_timing	
def validator_3(Q):
	# saveout = sys.stdout                                     
	# fsock = open('validator_3_match_data.csv', 'w')                             
	# sys.stdout = fsock
	queries = Q.keys()
	queries.sort(lambda x,y: cmp(Q[x][0].scan_number,Q[y][0].scan_number))
	scan_numbers = [Q[a][0].scan_number for a in queries]	
	D = zip(scan_numbers,queries)
	pairs = []
	isotopologue_groups = ["N","H","L","2H","2L"] # N = no isotopologue, H = heavy is an isotopologue, L = light is an isotopologue
	for i,line in enumerate(D):
		q1 = line[1] # query 
		s1 = line[0] # scan
		data1 = Q[q1][0] # all the mascot data
		low = max(s1-scan_width,0) # low scan number for test range
		while low not in scan_numbers and low <=s1: # find the first lower limit that is  real precursor scan number
			low += 1
		high = min(s1+scan_width, scan_numbers[-1])
		while high not in scan_numbers and high >= s1: # find the fist high scan number that is a precursor scan number
			high -= 1
		i1 = scan_numbers.index(low) 
		i2 = scan_numbers.index(high)
		R2 = D[i1:i2] # so this is the range of index elements to check
		if line in R2: R2.remove(line) # make sure to take out the one being tested
		for line2 in R2: # go through the range and start testing, compared to q1
			q2 = line2[1] # query
			s2 = line2[0] # scan
			data2 = Q[q2][0] # mascot data
			MW1 = data1.precursor_MW  # this is based on the measured precursor m/z and is = (precursor_mz * charge) - (charge * proton_mass)
			MW2 = data2.precursor_MW
			d = abs(MW1-MW2) # difference in mass (monoisotopic)
			# figure out which MW is the low and which is high
			MWmax = max(MW1,MW2)
			MWmin = min(MW1,MW2)
			# MW_test_group is setup in main and =  [['LYS', 6.0201260000000003], ['ARG', 10.00827], ['LYS_LYS', 12.040252000000001], ['ARG_ARG', 20.016539999999999], ['LYS_ARG', 16.028396000000001], ['LYS_ARG_ARG', 26.036666], ['ARG_ARG_LYS', 26.036666], ['LYS_LYS_LYS', 18.060378], ['ARG_ARG_ARG', 30.024809999999999]]
			''' This tests each scenario - each possible isotope (e.g. ARG, LYS, ARG-LYS) and one isotopologue on the light or heavy 
			    So, the three values in each line are N, L, H, and each line corresponds to a group in MW_test_group
			    Therefore the total number of tests for each pair of peptides will be 9*3 = 27 '''
			#MW_test = [   [            MWmax-a[1] - MWmin   ,    MWmax-a[1] - first_isotopologoue - MWmin ,   MWmax - a[1] - ( MWmin - first_isotopologue)           ] for a in MW_test_group]
			MW_test = [   [ MWmax-a[1] - MWmin,    
							MWmax-a[1] - first_isotopologue - MWmin,   
							MWmax - a[1] - ( MWmin - first_isotopologue), 
							MWmax-a[1] - second_isotopologue - MWmin, 
							MWmax - a[1] - ( MWmin - second_isotopologue)] for a in MW_test_group]

			PPM_test = [] # now for each MW found, figure out the PPM error
			for line in MW_test:
				newline = [a/MWmin*1000000 for a in line] # the ppm is calculated using the lower monoisotopic mass
				PPM_test.append(newline)
			# go through each one, see if it makes the grade and then figure out what the isotope and isotopologue were	
			for i,line in enumerate(PPM_test): # these are the isotopes
				for j,l in enumerate(line): # these are the isotopologues
					if abs(l)<ppm_cutoff:
						pme = l # remember this PME
						isotope = MW_test_group[i] # e.g. ARG, LYS, ARG-LYS, etc.
						isotopologue = isotopologue_groups[j] # H,L or N
						if MW1<MW2: # add on the new pais as light,heavy
							new_pair = [q1,q2]
						else:
							new_pair = [q2,q1]
						if new_pair not in [a[0] for a in pairs]:
							pairs.append([new_pair, isotope, isotopologue, abs(pme)])						
	#print "Q1, Q2, Scan 1, Scan 2, Mods 1, Mods 2, Peptide 1, Peptide 2, Score 1, Score 2, Isotope, Isotopologue, ABS PME, Match?"
	#pairs.sort(lambda a,b:cmp(abs(a[3]),abs(b[3])))	
	# for pair in pairs:
	# 	q1 = pair[0][0]
	# 	q2 = pair[0][1]
	# 	isotope = pair[1]
	# 	isotopologue = pair[2]
	# 	pme = pair[3]
	# 	s1 = Q[q1][0].scan_number
	# 	s2 = Q[q2][0].scan_number
	# 	line = [q1,q2,s1,s2,"'"+Q[q1][0].top_mods,"'"+Q[q2][0].top_mods,Q[q1][0].top_peptide,Q[q2][0].top_peptide, Q[q1][0].top_score,Q[q2][0].top_score,isotope[0], isotopologue, pme, Q[q1][0].top_peptide==Q[q2][0].top_peptide]
	# 	print ",".join([str(a) for a in line])
	# sys.stdout = saveout
	# fsock.close()
	
	qq = unique([a[0][0] for a in pairs]+[a[0][1] for a in pairs])
	pp = [p for p in qq if len(Q[p][0].peptides) > 0]
	if screen_print:
		print "Validator found",len(pairs)," unique pairs and",len(qq),"unique queries, of which",len(pp),"have peptides."
	return pairs

def corroborate_peptide(peptide,ions,mods):

 	sum_intensities = ions[3]
	# the ions coming in are for light or heavy. If iso_mass = 0, then we are light, if not 0 then we are heavy
	mod_string = ''
	
	for pos,mod in enumerate(mods):
		mod = int(mod)
		if mod <> 0:
			m = modifications[mod]
			key = m[0]
			mass = m[1]
			mod_string += str(mass)+'@'+str(pos+1)+" "

	f = PeptideFragment(1)
	f.analyze(peptide,1,mod_string)

	m = [a for a in f.peptide_mass()][1]
	
	b_ions_peptide = [a for a in f.b_ions()][:-1] + [m]
	y_ions_peptide = [m] + [a for a in f.y_ions()][1:]

	NH3_mass = 17.0306	
	b_ions_NH3_loss = []
	for a in range(1,len(peptide)+1):
		#print peptide,len(peptide),len(b_ions_peptide),a, iso_mass
		p = peptide[:a]
		m = b_ions_peptide[a-1]
		if 'R' in p or 'K' in p or 'Q' in p or 'N' in p:
			b_ions_NH3_loss.append(m-NH3_mass)
		else:
			b_ions_NH3_loss.append(0.0)
	y_ions_NH3_loss = []
	for a in range(len(peptide),0,-1):
		p = peptide[-a:]
		m = y_ions_peptide[len(peptide)-a]
		if 'R' in p or 'K' in p or 'Q' in p or 'N' in p:
			y_ions_NH3_loss.append(m-NH3_mass)
		else:
			y_ions_NH3_loss.append(0.0)

	H2O_mass = 18.0153
	b_ions_H2O_loss = []
	for a in range(1,len(peptide)+1):
		p = peptide[:a]
		m = b_ions_peptide[a-1]
		if 'S' in p or 'T' in p or 'E' in p or 'D' in p:
			b_ions_H2O_loss.append(m-H2O_mass)
		else:
			b_ions_H2O_loss.append(0.0)
	y_ions_H2O_loss = []
	for a in range(len(peptide),0,-1):
		p = peptide[-a:]
		m = y_ions_peptide[len(peptide)-a]
		if 'S' in p or 'T' in p or 'E'in p or 'D' in p:
			y_ions_H2O_loss.append(m-H2O_mass)
		else:
			y_ions_H2O_loss.append(0.0)

	proton_mass = 1.00727646688
	b_ions_pos_charge = []
	for a in range(1,len(peptide)+1):
		p = peptide[:a]
		m = b_ions_peptide[a-1]
		if 'R' in p or 'H' in p or 'V' in p:
			b_ions_pos_charge.append(  (m+ proton_mass ) / 2.0 )
		else:
			b_ions_pos_charge.append(0.0)
	y_ions_pos_charge = []
	for a in range(len(peptide),0,-1):
		p = peptide[-a:]
		m = y_ions_peptide[len(peptide)-a]
		if 'R' in p or 'H' in p or 'V'in p:
			y_ions_pos_charge.append( (m+ proton_mass ) / 2.0 )
		else:
			y_ions_pos_charge.append(0.0)

	all_ions = [b_ions_peptide, y_ions_peptide, b_ions_pos_charge, y_ions_pos_charge, b_ions_NH3_loss, y_ions_NH3_loss, b_ions_H2O_loss, y_ions_H2O_loss]	


	''' THIS PRINTS FRAGMENTATION TABLE '''

	if screen_print2:
		header = ["","b_ions",  "b_ions_pos_charge",  "b_ions_NH3_loss", "b_ions_H2O_loss", "y_ions", "y_ions_pos_charge", "y_ions_NH3_loss","y_ions_H2O_loss"]	
		D = zip(peptide,b_ions_peptide, b_ions_pos_charge, b_ions_NH3_loss,  b_ions_H2O_loss,y_ions_peptide, y_ions_pos_charge,  y_ions_NH3_loss,  y_ions_H2O_loss)
		table = [header] + D
		out = sys.stdout
		"PREDICTED FRAGMENTATION TABLE"
		pprint_table(out, table)
		print
		print
	

	score = [0,0]
	total_intensity = [0,0]
	for x,I in enumerate(all_ions):
		if x%2 == 0:
			z = 0 # B
		else:
			z = 1 # Y
		
		half_shift = False
		#if x ==2 or x ==3 :
		if x ==3 :
			half_shift = True
			
		#print

		for y,ion in enumerate(I):
			
			#print trunc(ion)
			# need to determine where we are in the peptide to know whether to look for non-shifters or shifters
			# for instance:
			#
			#      _____K______R
			#      _____K______R
			#
			#   In the first part, B will be non-shift, Y will shift, but in the second part, both will be shifted.
			p1 = peptide[:y]
			p2 = peptide[y:]
			#print p1,p2
			if ("K" not in p1 and "R" not in p1) and ("K" in p2 or "R" in p2):
				# B's will not shift
				# Y's will shift
				if z == 0:
					zz = 0 # look for B's in the non-shifter bin
				elif half_shift:
					zz = 2
				else:
					zz = 1 # look for Y's in the shifter bin
			elif ("K" in p1 or "R" in p1) and ("K" not in p2 and "R" not in p2): # like a c-terminal peptide
				# B's will shift
				# Y's will not shift
				if z == 0 and half_shift:
					zz = 2 # look for B's in the shifter bin
				elif z == 0:
					zz = 1
				else:
					zz = 0 # look for Y's in the non-shifter bin
			elif ("K" in p1 or "R" in p1) and ("K" in p2 or "R" in p2):
				# B's will shift
				# Y's will shift
				if half_shift:
					zz = 2
				else:
					zz = 1 # look in the shifter bin for everything
				
			#print "LEN",len(ions),zz
			pos = bisect.bisect_left([a[0] for a in ions[zz]],ion)
			#pos = bisect.bisect_left([ (a[0] +p_add) / half_shift_div for a in ions[zz]],ion)

			if pos>0 and ion > 0 and same_ion(ion,ions[zz][pos-1][0]):
				score[z] += 1
				total_intensity[z] += ions[zz][pos-1][1]
			elif pos<len(ions[zz]) and ion > 0 and same_ion(ion,ions[zz][pos][0]):
				score[z] += 1
				total_intensity[z] += ions[zz][pos][1]
			elif ion>0:
				all_ions[x][y] = 'X'
			else:
				all_ions[x][y] = '-'#ion * -1.0
								
	total_consecutive = 0

	for I in all_ions:
		#newmax = 1
		con = 0
		for i in range(0,len(I)-1):
			if I[i] <> "X" and I[i] <> "-" and I[i+1] <> "X" and I[i+1] <> "-" and i+1 <> len(I) - 1: # we have two consecutive in a row and we are not at the end
				#print I[i],I[i+1],i,len(I)
				con += 1
			elif I[i] <> "X" and I[i] <> "-" and I[i+1] <> "X" and I[i+1] <> "-": # we are at the end
				#print I[i],I[i+1],i,len(I)
				total_consecutive += con + 2
				con = 0
			elif I[i] <> "X" and I[i] <> "-" and (I[i+1] == "X" or I[i+1] == "-") and con>0: # we've reached the end of our run
				#print I[i], "END"
				total_consecutive += con + 1
				con = 0
	
	
	ions = zip(all_ions[0],all_ions[2],all_ions[4],all_ions[6],peptide,all_ions[1],all_ions[3],all_ions[5],all_ions[7])

	ions2 = zip(peptide,all_ions[0],all_ions[2],all_ions[4],all_ions[6],all_ions[1],all_ions[3],all_ions[5],all_ions[7])
	if screen_print2:
		table = [header] + ions2
		out = sys.stdout
		print "CORROBORATED FRAGMENTATION TABLE"
		pprint_table(out, table)
		print
		print "Total intensity |  total intensity / sum intensities"
		print total_intensity[0],total_intensity[0]/sum_intensities
		print total_intensity[1],total_intensity[1]/sum_intensities
		print
		print

	scoreb = sum([sum([(s>0 and s<>'X' and s<>'-') for s in a ]) for a in [b[0:4] for b in ions]] )
	scorey = sum([sum([(s>0 and s<>'X' and s<>'-') for s in a ]) for a in [b[5:9] for b in ions]] )
	
	percent_intensity = [total_intensity[0]/sum_intensities, total_intensity[1]/sum_intensities]
	

	return ions,total_consecutive,[scoreb,scorey],percent_intensity

def corroborate(D,Q):
	#total_consecutive = 0

	iso = D[2][0]
	peptides = [Q[a][0].top_peptide for a in D[0]]
	mods = [Q[a][0].top_mods[1:-1] for a in D[0]]
	output = []
	flag = 0
	if peptides[0] <> peptides[1]:
		flag = 1 # peptides are not the same
	'''D[1] are the ions. D[1][0] are the non-shift, D[1][1] are the shift, D[1][2] are the half-shifters, and D[1][3] are the non-matched and D[1][4] are the sum of ion intensities'''
	ionsL = [a[0] for a in D[1][0]],[a[0] for a in D[1][1]],[a[0] for a in D[1][2]],D[1][4][0]
	t1 = ionsL[0]
	t1 = unique(t1)
	t1.sort()
	t2 = ionsL[1]
	t2 = unique(t2)
	t2.sort()
	t3 = ionsL[2]
	t3 = unique(t3)
	t3.sort
	ionsL = (t1,t2,t3,ionsL[3])
	ionsH = [a[1] for a in D[1][0]],[a[1] for a in D[1][1]],[a[1] for a in D[1][2]],D[1][4][1]
	t1 = ionsH[0]
	t1 = unique(t1)
	t1.sort()
	t2 = ionsH[1]
	t2 = unique(t2)
	t2.sort()
	t3 = ionsH[2]
	t3 = unique(t3)
	t3.sort()
	ionsH = (t1,t2,t3,ionsH[3])


	ions = [ionsL, ionsH]
	peptideL = peptides[0]
	peptideH = peptides[1]
	peptides = [peptideL, peptideH]
	modsL = mods[0]
	modsH = mods[1]	
	mods = [modsL, modsH]
	# this checks to see if the predicticed isotope (eg. ARG, ARG-LYS) is consistent with the peptide, returns array of TRUE/FALSE, TRUE/FALSE
	# this works by checking to see if the expected number of K and then R are found in the peptide and returning an array for each. 
	# for example, if the isotope was LYS-ARG and the peptide was AAAAAKAAAR, it woudl return [TRUE, TRUE], but if the isotope was ARG, it would return [FALSE, TRUE]
	# this is the general case
	'''So, isotope check merely confirms that the isotope used for the mass difference is appropriate for the peptide guess.
	   It remins to be seen whether the mods are consistent with a heavy or light - still need to check this.'''
	isotope_check = [str(a) for a in   [sum([[mod_names[[b[0] for b in mod_names].index(a)][1] for a in D[2][0].split('_')].count(m) == peptide.count(m) for m in [a[1] for a in mod_names]])==len(mod_names) for peptide in [peptideL,peptideH]]]
	
	mods_check = [0,0]
	for i,line in enumerate(zip([modsL,modsH],[peptideL,peptideH])):
		# line is mods, peptide. i indicates light (0) or heavy (1)
		mods_ = line[0]
		peptide = line[1]
		if not isotope_check[i]:
			mods_check[i] == 'False'
		else:
			mod_indices = [peptide.index(a) for a in peptide if a in [b[1] for b in mod_names]]
			s = sum([int(mods_[a]) for a in mod_indices]) 
			if (s <> 0) == i:
				mods_check[i] = 'True'
			else:
				mods_check[i] = 'False'
		
	output = []

	def rev_mods(m,p):
		''' reverses mods.If K or R, 0 becomes 4 or 5 and 4 or 5 becomes 0'''
		new = ""
		for line in zip(m,p):
			aa = line[1]
			mod = line[0]
			if aa not in [a[1] for a in mod_names]:
				new += mod
			else:
				if mod == str(mod_dict[aa]):
					new += '0'
				else:
					new += str(mod_dict[aa])
		return new

	if peptideL == peptideH and peptideL <> "None": # this is a match, so display the fragmentation patterns and count the number of matches
		output.append(      [  corroborate_peptide(peptideL,ions[0],mods[0]), corroborate_peptide(peptideL,ions[1],mods[1])  ]  )
	elif (isotope_check[0] == 'True' and mods_check[0] == 'True') or (isotope_check[1] == 'True' and mods_check[1] == 'True'): 
		if (isotope_check[0] == 'True' and mods_check[0] == 'True'):# the mods seem to work, so the other peptide might actually be an incorrect assignment	
			new_mods = rev_mods(mods[0],peptideL)
			output.append(  [  corroborate_peptide(peptideL,ions[0],mods[0]),corroborate_peptide(peptideL,ions[1],new_mods) ] ) # need to fix the second mods
		if (isotope_check[1] == 'True' and mods_check[1] == 'True'): 
			new_mods = rev_mods(mods[1],peptideH)
			output.append(  [  corroborate_peptide(peptideH,ions[0],new_mods),corroborate_peptide(peptideH,ions[1],mods[1]) ] ) # need to fix the first mods
	return output,isotope_check,mods_check

def same_ion_size_and_intensity(ion1,ion2,max_intensity):
	'''ion1 and ion2 are [m/z,intensity] and max is max intensity for each [max,max] and threshold is min_ion_percent_intensity and max_ion_difference (globals)'''
	m1 = ion1[0]
	m2 = ion2[0]
	i1 = ion1[1]
	i2 = ion2[1]
	max1 = max_intensity[0]
	max2 = max_intensity[1]
	ppm_calc = calc_ppm(m1,m2)
	'''this line ensures that the two m/z are the same, that the peaks heights are above minimum (eg. 1% of max) and that the peak heigh difference isn't too much '''
	rel_int_1 = i1/max1
	rel_int_2 = i2/max2
	if (abs(ppm_calc) < ms2_ppm_tolerance) and (abs(rel_int_1 - rel_int_2) / max(rel_int_1, rel_int_2)) < max_ion_difference:
		return True
	else:
		return False

def validator_3_function(pairs,Q):
	
	'''
	This function will take each pair and bin the MS2 ions as B and Y. 
	Formerly, this was done by looking for shifters and non-shifters.
	This function was modified to now account for missed cleavages - so the B and Y ions might have one (at least one for Y) internal isotopes
	'''
	all_matches = []
	for pair in pairs:

		q1 = pair[0][0]
		q2 = pair[0][1]
		c = [[],[],[]] # non-shifter, shifters and half-shifter matches
		iso = pair[1][0]
		d1 = Q[q1][0]
		d2 = Q[q2][0]
		MS2_1 = d1.MS2_ions
		MS2_2 = d2.MS2_ions
		MS2_ions_1 = [a for a in MS2_1] # just grab the m/z
		MS2_ions_2 = [a for a in MS2_2]		
		MS2_ions_1.sort(lambda a,b:cmp(a[0],b[0]))
		MS2_ions_2.sort(lambda a,b:cmp(a[0],b[0]))
		max_intensity_1 = max([a[1] for a in MS2_ions_1])
		max_intensity_2 = max([a[1] for a in MS2_ions_2])
		max_intensities = [max_intensity_1, max_intensity_2]
		
		if screen_print:
			print "MS1 cutoff", MS1_cutoff
			print "MS2 cutoff", MS2_cutoff
			print "QUERIES",q1,q2
			print "SCANS",d1.scan_number, d2.scan_number
		
			s = ["LIGHT","HEAVY"]
			print "These are all the MS2 ions"
			for x,ions in enumerate([MS2_ions_1, MS2_ions_2]):
				print
				print s[x]
				print '%8s'%("m/z"),'%8s'%("intensity") ,'%8s'%("percent max") 
				for ion in ions:
					print '%8.3f'%(ion[0]),'%8.2f'%(ion[1]), '%8.4f'%(ion[1]/max_intensities[x])
				print
				print
				print
		
		#
		#  CONSIDER MODIFYING THIS SECTION TO CHANGE ALGORITHM ON FILTERING OUT LOW IONS
		#
		'''Filter out the ions that do not make the intensity criteria
		   THIS IS AN IMPORTANT STEP. A lot of ions will be lost here.'''
		MS2_ions_1 = [a for a in MS2_ions_1 if a[1]/max_intensity_1>min_ion_percent_intensity]
		MS2_ions_2 = [a for a in MS2_ions_2 if a[1]/max_intensity_2>min_ion_percent_intensity]
		MS2 = [a[0] for a in MS2_ions_2]

		if screen_print:
			#
			#  PRINTOUT MS2 ions
			#
			print "A filter was applied to remove those ions not greater than",min_ion_percent_intensity,"of the max intensity"
			print "The max intensity for light and heavy are"
			print max_intensity_1, min_ion_percent_intensity * max_intensity_1
			print max_intensity_2, min_ion_percent_intensity * max_intensity_2
		
			print "These are the ions that are left after filtering."
			for x,ions in enumerate([MS2_ions_1, MS2_ions_2]):
				print "These are all the MS2 ions"
			
				print '%8s'%("m/z"),'%8s'%("intensity") ,'%8s'%("percent max") 
				for ion in ions:
					print '%8.3f'%(ion[0]),'%8.2f'%(ion[1]), '%8.4f'%(ion[1]/max_intensities[x])
				print
				print
				print

		M = [max_intensity_1, max_intensity_2]
		S = [ sum(a[1] for a in MS2_ions_1), sum(a[1] for a in MS2_ions_2)]
		isotopologue = pair[2]
		PME = pair[3]
		''' Look for shifters and non-shifters
			For internal missed cleavages, will need to consider that B ions might be shifting and Y ions might be single or double (or more) shifted '''		

		isotopes = iso.split('_')			
		
		if screen_print:
			print "These are the isotopes",isotopes
			print "And these are the masses we will be checking."
			for n in range(0,len(isotopes)+1): # need to consider the scenarios when the difference between ions in 0, 1... up to the full number of isotopes - 
				combos = itertools.combinations(isotopes,n) # get the possible cobinations from the isotopes  - so if ARG_LYS_ARG, gets ARG, LYS, ARG_LYS, etc.
				combos = [a for a in combos]
				combos = unique(combos) 

				for combo in combos:
					m = 0 # start with adding a mass of 0. the first combo will be the empty set anyway. This will find the b ions for the non-shifters
					for isotope in combo: # sequentially add each isotope mass
						m += isotope_masses[isotope]
					print "Mass to check",m
			print
			print


		''' Group the ions into shifters and non-shifters '''
		for ion in MS2_ions_1: # consider each MS2 in sequence from the first group
			T = len(MS2_ions_2)
			for n in range(0,len(isotopes)+1): # need to consider the scenarios when the difference between ions in 0, 1... up to the full number of isotopes - 
				combos = itertools.combinations(isotopes,n) # get the possible combinations from the isotopes  - so if ARG_LYS_ARG, gets ARG, LYS, ARG_LYS, etc.
				combos = [a for a in combos]
				combos = unique(combos) 
			
				for combo in combos:
					# if len(combo) == 0: # non-shifter
					# 	x = 0
					# else:
					# 	x = 1
					
					m = 0.0 # start with adding a mass of 0. the first combo will be the empty set anyway. This will find the b ions for the non-shifters
					for isotope in combo: # sequentially add each isotope mass
						m += isotope_masses[isotope]
					if m == 0.0: # non-shifter
						X =  [ [0,0.0] ]# looking for non-shifters. First number is position in c (non-shift, shift, half-shift)
					else:
						X  = [ [1,m],[2,m/2.0]     ] # this sets up to look for shifters and half shifters and add to the c array accordingly
						
					for XX in X: # for non-shifters, this is only one item- the x will be 0 and the m will be 0. But for shifters, this will be two items - first for full shift and second for half shift
						x = XX[0] # place in array (0 = non-shift, 1 = shift, 2 = half-shift)
						m = XX[1] # the shifted amount to look for 
						
						pos = bisect.bisect_left(MS2,ion[0]+m) # find the point of insertion in the group of MS2 ions from the heavy
						#print "XXX",x,m
						if pos>0 and same_ion_size_and_intensity([ion[0]+m,ion[1]],MS2_ions_2[pos-1],M):
							c[x].append([ion,MS2_ions_2[pos-1]]) # tack the match onto the end of B or Y or both						
						elif pos<T and same_ion_size_and_intensity([ion[0]+m,ion[1]],MS2_ions_2[pos],M):
							c[x].append([ion,MS2_ions_2[pos]])
					
						
		not_matched_light = [a for a in MS2_ions_1 if (a not in [b[0] for b in c[0]] and a not in [b[0] for b in c[1]]) and a[1]/max_intensity_1 > min_ion_percent_intensity]
		not_matched_heavy = [a for a in MS2_ions_2 if( a not in [b[0] for b in c[0]] and a not in [b[0] for b in c[1]]) and a[1]/max_intensity_2 > min_ion_percent_intensity]
		if screen_print:
			print "Max Intensities"
			print '%8.3f'%(M[0]),	 '%8.3f'%(M[1])
			header = ["m/z light","m/z heavy","diff m/z","intensity L","intensity H","normalized intensity L","normalized intensity H","Diff normalized intensity","Relative Intensity difference"]
			out = sys.stdout
			print "Non-shifters"
			d = []
			for line in c[0]:
				d.append(['%8.3f'%(line[0][0]),  '%8.3f'%(line[1][0]), '%8.3f'%(abs(line[1][0]-line[0][0])), '%8.3f'%(line[0][1]), '%8.3f'%(line[1][1]), '%8.3f'%(line[0][1]/M[0]), '%8.3f'%(line[1][1]/M[1]),  '%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])),'%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])  / max( [line[0][1]/M[0],line[1][1]/M[1]]  ))]) 
			table = [header] + d
			pprint_table(out, table)
			print
			print "Shifters"
			d = []
			for line in c[1]:
				d.append(['%8.3f'%(line[0][0]),  '%8.3f'%(line[1][0]), '%8.3f'%(abs(line[1][0]-line[0][0])), '%8.3f'%(line[0][1]), '%8.3f'%(line[1][1]), '%8.3f'%(line[0][1]/M[0]), '%8.3f'%(line[1][1]/M[1]),  '%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])),'%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])  / max( [line[0][1]/M[0],line[1][1]/M[1]]  ))]) 
			table = [header] + d
			pprint_table(out, table)
			print
			print "Half-Shifters"
			d = []
			for line in c[2]:
				d.append(['%8.3f'%(line[0][0]),  '%8.3f'%(line[1][0]), '%8.3f'%(abs(line[1][0]-line[0][0])), '%8.3f'%(line[0][1]), '%8.3f'%(line[1][1]), '%8.3f'%(line[0][1]/M[0]), '%8.3f'%(line[1][1]/M[1]),  '%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])),'%8.3f'%(abs(line[0][1]/M[0]-line[1][1]/M[1])  / max( [line[0][1]/M[0],line[1][1]/M[1]]  ))]) 
			table = [header] + d
			pprint_table(out, table)
			print
		
		c.append([not_matched_light,not_matched_heavy]) # c is [[non-shiftersL, non-shiftersH], [shiftersL, shiftersH], [half-shiftersL, half-shiftersH], [not-matchedL, not-matchedH], [sumions1,sumions2]]		
		c.append(S)		
		mass_1 = d1.precursor_MW
		mass_2 = d2.precursor_MW
		final_score = 100.0  * float(len(c[0]) + 2 * len(c[1])) / mass_1 # final score uses number of non-shifters and shifters
		all_matches.append([[q1,q2],c,[iso,isotopologue,PME], final_score])
	return all_matches

#@print_timing
def validator_studies():

	@print_timing
	def validator_3_analysis():
		if screen_print:
			print
			print "##################################"
			print "#########   VALIDATOR 3 ##########"
			print "##################################"
			print
		pairs_3 = validator_3(query_dictionary) # [['7733', '7753'], ['ARG', 10.00827], 'N', -0.12918151875875858]
		if screen_print:
			print
			print "##################################"
			print "#########   VALIDATOR 3e ##########"
			print "##################################"
			print		
		matches_3 = validator_3_function(pairs_3,query_dictionary) # get scores, each line is ([[q1,q2],c,[iso,isotopologue,PME], final_score]) where c has the b ions, y ions and non-matches

		matches_3.sort(lambda a,b:cmp(b[3],a[3])) # sort by final score (B + 2Y for instance)
		Q = query_dictionary
		peptides_matched_light_heavy = []
		proteins_matched_light_heavy = []
		peptides_mismatched_light_heavy = []
		proteins_mismatched_light_heavy = []
		proteins_matched_any = 0
		
		# these are the column headings and column widths. The second header is the name of the category for the fragmentation sheet
		output = [[	["Pair",7, 0, "Pair" ],
					["Score [(S+2NS)/Mass]",12, 0, "Score [(S+2NS)/Mass]" ],
					["Query L",10, 0, "Query"],
					["Query H",10, 1],
					["Scan L",10, 0, "Scan"],
					["Scan H",10, 1],
					['Precursor m/z L',15,0, "Precursor m/z"],
					['Precursor m/z H',15,1],
					['Charge L',9,0,"Charge"],
					['Charge H',9,1],
					['Charges match',10,0, "Charges match?"],
					['Neutral mass L',15,0, "Neutral mass"],
					['Neutral mass H',15,1],
					["Isotope",15, 0, "Isotope"],
					["Isotopologue",15, 0, "Isotopologue"],
					["ABS PME",10, 0, "ABS PME"],
					["Non-Shifters",10, 0, "Non-Shifters / Shifters"],
					["Shifters",10, 1],
					["Half-Shifters",10, 1],
					["Unmatched L",10, 0, "Unmatched L/H"],
					["Unmatched H",10, 1],
					["Peptide L",35, 0, "Peptide"],
					["Peptide H",35, 1],
					["Isotope match L",12, 0, "Isotope match?"],
					["Isotope match H",12, 1],
					["Pep found L",15, 0, "Peptide found"],
					["Pep found H",15, 1],					
					['Peptide calc mz L', 18,0,"Peptide calc m/z"],
					['Peptide calc mz H', 18,1],
					["Peptide PME L",15, 0, "Peptide PME"],
					["Peptide PME H",15, 1],					
					["Match?",7, 0, "Match?"],
					["Mods L",35, 0, "Mods"],
					["Mods H",35, 1],
					["Mods match L",12, 0, "Mods match?"],
					["Mods match H",12, 1],
					["Isotope + Mods match L",12, 0, "Isotope + Mods match?"],
					["Isotope + Mods match H",12, 1],
					["Mascot score L",14, 0, "Mascot score"],
					["Mascot score H",14, 1],
					["Mascot proteins L",22, 0, "Protein"],
					["Mascot proteins H",22, 1],
					["Total ions L",15,0,"Total ions"],
					["Total ions H",15,1],
					["Bs L",10, 0, "Ions matched light"],
					["Ys L",10, 1],
					["Sum ions L",13, 0, "Sum light"],
					["Consecutive L",15, 0, "Consecutive"],
					["B L % coverage",20, 0, "B/Y L % coverage"], 
					["Y L % coverage",20, 1], 
					["B+Y L % coverage",20, 0,"B+Y L % coverage"],
					["Bs H",10, 0, "Ions matched heavy"],
					["Ys H",10, 1],
					["Sum ions H",13, 0, "Sum heavy"],
					["Consecutive H",15, 0, "Consecutive"],
					["B H % coverage",20, 0, "B/Y H % coverage"], 
					["Y H % coverage",20, 1], 
					["B+Y H % coverage",20, 0,"B+Y H % coverage"],
					["B L/H Matches",20, 0, "B/Y L/H Matches"], 
					["Y L/H Matches",20, 1], 
					["B/Y Match Score",20, 0,"B/Y Match score"],
					["Bs L Alt",10, 0, "Ions matched light alt"],
					["Ys L Alt",10, 1],
					["Sum ions L Alt",18, 0, "Sum light alt"],
					["Consecutive L Alt",20, 0, "Consecutive alt"],
					["B L alt % coverage",20, 0, "B/Y L alt % coverage"], 
					["Y L alt % coverage",20, 1], 
					["B+Y L alt % coverage",20, 0,"B+Y L alt % coverage"],
					["B L/H Matches L L",20, 0, "B/Y L/H Matches L L"], 
					["Y L/H Matches L L",20, 1], 
					["B/Y Match Score L L",20, 0,"B/Y Match score L L"], 
					["Bs H Alt",10, 0, "Ions matched heavy alt"],
					["Ys H Alt",10, 1],
					["Sum ions H Alt",18, 0, "Sum heavy alt"],
					["Consecutive H Alt",20, 0, "Consecutive alt"],
					["B H alt % coverage",20, 0, "B/Y H alt % coverage"], 
					["Y H alt % coverage",20, 1], 
					["B+Y H alt % coverage",20, 0,"B+Y H alt % coverage"],
					["B L/H Matches H H",20, 0, "B/Y L/H Matches H H"], 
					["Y L/H Matches H H",20, 1], 
					["B/Y Match Score H H",20, 0,"B/Y Match score H H"] ]]
					
		# these are peptides that are found in pairs
		good_peptides_found = [Q[a[0][0]][0].top_peptide for a in matches_3 if Q[a[0][0]][0].top_peptide == Q[a[0][1]][0].top_peptide and  Q[a[0][0]][0].top_peptide <> "None"] 
		good_peptides_found = unique(good_peptides_found)
		good_peptides_found.sort()
		
		for X,line in enumerate(matches_3):
			q1 = line[0][0]
			q2 = line[0][1]
			D1 = Q[q1][0]
			D2 = Q[q2][0]

			pep = [D1.top_peptide,D2.top_peptide]
			pep.sort()
			pro = [D1.top_protein,D2.top_protein]
			pro.sort()
			
			#################################################
			#												#
			#    COUNT UP PEPTIDES AND PROTEINS MATCHED   	#
			#												#
			#################################################
			
			if pep[0] == pep[1]:
				peptides_matched_light_heavy.append(pep)
			else:
				peptides_mismatched_light_heavy.append(pep)

			if pro[0] == pro[1]:
				proteins_matched_light_heavy.append(pro)
			else:
				proteins_mismatched_light_heavy.append(pro)
			
			flag = 0
			for protein in pro[0]:
				if protein in pro[1]:
					proteins_matched_any += 1
					flag = 1
					break
			if flag == 0:
				for protein in pro[1]:
					if protein in pro[0]:
						proteins_matched_any += 1
						break

			
			#################################################
			#												#
			#   			 GENERATE OUTPUT 			  	#
			#												#
			#################################################
			
			'''Checks to see if the mods listed are correct for the peptide.'''
			#isotope_check = [str(a) for a in   [sum([[mod_names[[b[0] for b in mod_names].index(a)][1] for a in line[2][0].split('_')].count(m) == peptide.count(m) for m in [a[1] for a in mod_names]])==len(mod_names) for peptide in [D1.top_peptide,D2.top_peptide]]]
			ions_set, isotope_check, mods_match = corroborate(line,Q) # need to figure out if either of these two peptides matches better
			# ions_set[0] = ions, ions_set[1] = consecutive, ions_set[2] = scoreb,scorey, ions_set[3] = intensity_coverage[b],intensity_coverage[y]
			
			isotope_mods_match = [0,0]
			# print "I",isotope_check
			# print "M",mods_match
			for i,l in enumerate(zip(isotope_check,mods_match)):
				isotope_mods_match[i] = str((l[0] == 'True' and l[1] == 'True'))
			isotope_and_mods_match_L = isotope_mods_match[0]
			isotope_and_mods_match_H = isotope_mods_match[1]
			# print isotope_mods_match
			# print
			# print
			pair =  X+1
			score = line[3]
			query_light = int(line[0][0])
			query_heavy = int(line[0][1])
			scan_light =  D1.scan_number
			scan_heavy =  D2.scan_number
			precursor_L = D1.precursor_mz
			precursor_H = D2.precursor_mz
			charge_light = D1.charge
			charge_heavy = D2.charge
			charges_match = str(charge_light == charge_heavy)
			neutral_mass_light = D1.precursor_MW
			neutral_mass_heavy = D2.precursor_MW
			isotope =  line[2][0]
			isotopologue = line[2][1]
			abs_PME =  line[2][2]
			non_shifters =  len(line[1][0])
			shifters = len(line[1][1])
			half_shifters = len(line[1][2])
			unmatched_light = len(line[1][3][0])
			unmatched_heavy = len(line[1][3][1])
			peptide_light = D1.top_peptide
			peptide_heavy = D2.top_peptide
			if D1.top_peptide <> D2.top_peptide:				
				peptide_found_light = str(D1.top_peptide in good_peptides_found)
				peptide_found_heavy = str(D2.top_peptide in good_peptides_found)
			else:
				peptide_found_light = ""
				peptide_found_heavy = ""				
			peptide_mz_light = D1.top_calc_mz
			peptide_mz_heavy = D2.top_calc_mz
			peptide_PME_L = D1.top_PME
			peptide_PME_H = D2.top_PME
			peptide_match = str((peptide_light == peptide_heavy and peptide_light <> 'None'))
			mods_light =  D1.top_mods[1:-1]
			mods_heavy =  D2.top_mods[1:-1]
			isotope_check_light  = isotope_check[0]
			isotope_check_heavy = isotope_check[1]
			mods_match_light  = mods_match[0] # mods match means the light or heavy mods are correct (eg. 0 for light, not 0 for heavy)
			mods_match_heavy = mods_match[1]
			mascot_score_light =  D1.top_score
			mascot_score_heavy =  D2.top_score
			if D1.top_protein <> "None":
				mascot_proteins_light =  " / ".join(D1.top_protein)
			else:
				mascot_proteins_light = "-"
			if D2.top_protein <> "None":
				mascot_proteins_heavy = " / ".join(D2.top_protein)
			else:
				mascot_proteins_heavy = "-"
			total_ions_L = len(D1.MS2_ions)
			total_ions_H = len(D2.MS2_ions)
	
			
			NL = [	pair, score, query_light, query_heavy, scan_light, scan_heavy, precursor_L, precursor_H, charge_light, charge_heavy, charges_match, neutral_mass_light, neutral_mass_heavy,
					isotope, isotopologue, abs_PME, non_shifters, shifters, half_shifters,
					unmatched_light, unmatched_heavy, peptide_light, peptide_heavy, isotope_check_light, isotope_check_heavy,
					peptide_found_light, peptide_found_heavy, peptide_mz_light, peptide_mz_heavy, peptide_PME_L, peptide_PME_H, peptide_match,
					mods_light, mods_heavy, mods_match_light, mods_match_heavy, isotope_and_mods_match_L, isotope_and_mods_match_H,
					mascot_score_light, mascot_score_heavy, mascot_proteins_light, mascot_proteins_heavy, total_ions_L, total_ions_H]
					
			output.append([NL,ions_set])

		writeexcel(output,output_xls)		
		queries_target = unique([a[0] for a in [b[0] for b in matches_3]] + [a[1] for a in [b[0] for b in matches_3]])
		pairs_target = unique([a[0] for a in matches_3])
		
		# generate information file
		saveout = sys.stdout                                     
		#fsock = open('information_file.txt', 'w')
		fsock = open(output_txt,'w')                             
		sys.stdout = fsock
		print "ppm cutoff:",ppm_cutoff
		print "scan width:",scan_width
		print "ms2_ppm_tolerance:",ms2_ppm_tolerance
		print "min_ion_percent_intensity:",min_ion_percent_intensity
		print "max_ion_difference:",max_ion_difference
		print
		print "Number of unique queries:", len(queries_target)
		print "Number of unique pairs:", len(pairs_target)
		print "Number of pairs where peptides are same:",  len(unique([a for a in pairs_target if Q[a[0]][0].top_peptide == Q[a[1]][0].top_peptide]))
		print "Number pairs with matched peptides:",len(peptides_matched_light_heavy),"(Unique = "+str(len(unique(peptides_matched_light_heavy)))+")"
		print "Number pairs with mismatched peptides:",len(peptides_mismatched_light_heavy),"(Unique = "+str(len(unique(peptides_mismatched_light_heavy)))+")"
		print "Number pairs with matched proteins:",len(proteins_matched_light_heavy),"(Unique = "+str(len(unique(proteins_matched_light_heavy)))+")"
		print "Number pairs with mismatched proteins:",len(proteins_mismatched_light_heavy),"(Unique = "+str(len(unique(proteins_mismatched_light_heavy)))+")"
		print "Number of queries where any proteins matched:",proteins_matched_any
		sys.stdout = saveout                                     
		fsock.close()                                            
		
		queries = [queries_target]
		t_v3e = unique(queries[0])
		print_above_and_below_35(t_v3e,query_dictionary)
		print 'Done with analysis - check your output file.'
		return
	validator_1(query_dictionary)	
	validator_3_analysis()
	return

#@print_timing
def main(pickled = False):
	global query_dictionary
	global modifications
	global mods_to_check
	global poss_mods
	global MW_test_group
	global mod_names
	global isotope_masses
	global mod_dict
	current_dir = os.getcwd()

	filenames = [input_file]
	for filename in filenames:
		print "Starting",filename
		#os.chdir(current_dir)
		data = readfile(filename)
		dirname = os.path.dirname(filename)
		if not os.path.exists(dirname):
			os.makedirs(dirname)
		new_data = []
		# load in the whole .DAT file
		for line in data:
			if '\r' in line:
				line = line[:-1]
			new_data.append(line)
		data = new_data
		if screen_print:
			print "Data file has",len(data),"lines."
		content = make_data_dictionary(data)
		sections = content.keys()		
		mass_data = content['masses']
		modifications = {0:['N',0.0]} # instantiate the dictionary of modifications
		# find the modification information in the DAT file
		for line in mass_data:
			if 'delta' in line:
				n = int(line.split('=')[0].split('delta')[1])
				m = float(line.split('=')[1].split(',')[0])
				aa = line.split(' ')[1].split('(')[1].split(')')[0]
				if aa == "N-term":
					aa = 1
				elif aa == "K":
					lys_mod = n
					delta_lys6 = m
				elif aa == "R":
					arg_mod = n
					delta_arg10 = m
				modifications[n] = [aa,m]
		# will need to generalize this part for other isotope labeling
		# this script is currently only written for LYS/ARG SILAC
		if not('lys_mod' in dir() or 'arg_mod' in dir()):
			print filename, "not SILAC"
			continue # with next file in list
		mods_to_check = [lys_mod, arg_mod]
		poss_mods = [delta_lys6, delta_arg10]
		MW_test_group = []
		mod_names = [["LYS",'K'],["ARG","R"]]
		mod_dict = {} # will print corresponding number for amino acid
		for x,m in enumerate([a[1] for a in mod_names]):
			mod_dict[m] = mods_to_check[x]
		for m in "ABCDEFGHIJLMNOPQSTUVWXYZ":
			mod_dict[m] = 0
		# possible combinations of mods allowed (e.g. ARG-LYS-LYS)
		combos = ['0','1','00','11','01']#,'011','110','000','111', '0111', '1000', '1100', '1111', '0000']
		for c in combos:
			MW_test_group.append(["_".join([mod_names[int(b)][0] for b in c]),sum([poss_mods[int(cc)] for cc in c])])
		isotope_masses = {}
		# pre-calculate all the possible isotope masses
		for line in MW_test_group:
			key = line[0]
			value = line[1]
			isotope_masses[key] = value
		queries = [content[a] for a in content if 'query' in a]

		if not pickled: # if we are not using a smaller, pickled testing data set
			query_dictionary = make_query_dictionary(content["summary"],content["peptides"],queries)
			if screen_print:
				print "...done making target dictionary"			
			keys = query_dictionary.keys()
			# uncomment these lines in order to make a smaller testing set and then set the flag in the definition line to pickled = True
		 	#query_dictionary = extract(query_dictionary,[input_query_1,input_query_2])
		
			#query_dictionary_small = extract(query_dictionary,keys[3000:4000])	
			# query_dictionary_small = extract(query_dictionary,[input_query_1,input_query_2])		
			# pickle_data(query_dictionary_small,"query_dictionary_small.pickle")
			# exit()
			#pickle_data(query_dictionary_small,"query_dictionary_small.pickle")
			#exit()
		else:
			print current_dir + '/'  + "query_dictionary_small.pickle"
			query_dictionary = load_pickle_data(current_dir + '/'  + "query_dictionary_small.pickle")
			#query_dictionary = load_pickle_data("query_dictionary_small.pickle")
		validator_studies() # runs Validator 1 and 3 and 3e in succession0
	return

main()
subprocess.call(["perl","ValidatorParse_cli.pl",output_folder, output_folder+"/"+xls_name])
if not cpas_file == '':
	#pairs.csv CrossRef and Quant
	subprocess.call(["perl","ValidatorCrossRef.pl",output_csv, cpas_file, output_folder+"/crossref_pass.txt", "Val1"])
	subprocess.call(["perl", "HanashQuantRevise.pl", output_folder+"/crossref_pass.txt", "all"])
	
	for parse_type in ["MascotMatch", "HLRecover", "AllMatch"]:
		parse_output = output_folder+"/"+parse_type+"/"
		subprocess.call(["perl","ValidatorCrossRef.pl",parse_output+xls_name.split(".")[0]+"_"+parse_type+".tsv", cpas_file, parse_output+"crossref_pass.txt", "Val3"])
		subprocess.call(["perl", "HanashQuantRevise.pl", parse_output+"crossref_pass.txt", "all"])
	

print "FINISHED"

