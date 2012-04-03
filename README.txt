The following describe the scripts used to generate the data for our manuscript. The details for analyzing one slice are below. Sample data are included, along with the respective output files.

VALIDATOR

Running Validator requires several prerequisites.
A minimum of Python 2.6 is required. Earlier versions may work, but the software has only been tested using Python 2.6 and 2.7.

A C library has been written to speed up calculations of peptide fragments. The code is located in code/fragment-singleton and must be compiled separately for each architecture and python version. A sample makefile is included. Enter the path to your python include folder in the first line. You can sometimes find this by using the command:
python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
or by issuing python -v from the command line.
The two resulting files, "libpeptidefragmentsingleton.so" and "PeptideFragmentSingleton.so" must be in the same directory as the validatory.py script.

Other requirements are the following modules: cPickle, itertools, xlwt. 
Padnums.py is also required if using the screen-print function. Uncomment the import line to use this functionality. Padnums.py must be in the same folder as the validator.py script.

To run validator, first copy or move the folder containing your Mascot DAT files to the data folder. Then change to the code directory and execute the script as follows:

$ python validator.py <data folder name>

e.g. $ python validator.py "sample_data"

The program will take about five minutes to run for a 100mb DAT file.

The resulting files will be placed in "output_files/" in a folder named the same as your data folder. There are three data files produced:

information_file.txt contains some information about the variables used in the run
Val_1_pairs.csv contains the list of Validator 1 pairs
validator_3.xls contains the validator 3 data

If desired, one could stop here and study the output. 

CROSSREF

The inputs to this script are (1) Val_1_pairs.csv and (2) a file exported from CPAS/Peptide prophet. Export with no filters and heavy area > 0 and light area > 0. 

To run the crossref program, change to the code directory execute the command using the following format:

$ perl ValidatorCrossRef.pl <val_1_file> <CPAS file> <output_filename> Val1

e.g. perl ValidatorCrossRef.pl ../output_files/sample_data/F004642_sec_12/Val_1_pairs.csv ../supporting_data/MS2Runs_E1sec12.txt ../output_files/sample_data/F004642_sec_12/E1sec12_crossref_pass.txt Val1

You should now have a file in output_files/sample_data/<your folder name> with the name you supplied above.


QUANT

To run the quant script, you must have the t-table and p-table files located in supporting_data. These have been supplied for you.
The syntax of the program is as follows:

$ perl HanashQuantRevise.pl <name of cross ref file> all

e.g. perl HanashQuantRevise.pl ../output_files/sample_data/F004642_sec_12/E1sec12_crossref_pass.txt all

As before, you must do this from within the code folder. You will now have two additional files in your output folder - all.events and all.quant.

VALIDATOR PARSE

The script ValidatorParse.pl has an important prerequisite, "Spreadsheet-ParseExcel" at least version 0.59." This module has several pre-requisites of its own. These files are located in the folder code/perl_modules_required.

Before running the script, create a file called ValidatorList in the supporting_files folder that contains the path to excel output file from Validator. This path should be relative to the data folder. On the second line, list he name of the excel output file.

The syntax for the script is:

$ perl ValidatorParse.pl ../supporting_data/ValidatorList

Three files will be created:
AllMatch
HLRecover
MascotMatch




 