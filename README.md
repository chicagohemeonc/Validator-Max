# VALIDATOR
### Requirements
* Python2.6+
* cPickle
* itertools
* xlwt
* pybindgen

### Installation

```
git submodule update --init
cd code/fragment-singleton
python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
```
Copy the output to the PYTHON_INCLUDE in Makefile
```
make
cp *.so ../.
```
Set environment variable LD_LIBRARY_PATH="."

### USAGE - Piecemeal
All code is run from code/
#### Validator
* Copy or move the folder containing your Mascot DAT files to the data folder

```
python validator.py <data folder name>
# eg. python validator.py "sample_data"
```

The program will take about five minutes to run for a 100mb DAT file.

The resulting files will be placed in "output_files/" in a folder named the same as your data folder. There are three data files produced:
* information_file.txt contains some information about the variables used in the run
* Val_1_pairs.csv contains the list of Validator 1 pairs
* validator_3.xls contains the validator 3 data

#### CROSSREF

The inputs files to this script are:
* Val_1_pairs.csv
* A file exported from CPAS/Peptide prophet. Export with no filters and heavy area > 0 and light area > 0. 

```
perl ValidatorCrossRef.pl <val_1_file> <CPAS file> <output_filename> Val1
# e.g. perl ValidatorCrossRef.pl ../output_files/sample_data/F004642_sec_12/Val_1_pairs.csv ../supporting_data/MS2Runs_E1sec12.txt ../output_files/sample_data/F004642_sec_12/E1sec12_crossref_pass.txt Val1
```

The resulting files will be placed in <output_filename> with the name you supplied above.

#### QUANT
The inputs files to this script are t-table and p-table files located in supporting_data. These have been supplied for you.

```
perl HanashQuantRevise.pl <name of cross ref file> all
# e.g. perl HanashQuantRevise.pl ../output_files/sample_data/F004642_sec_12/E1sec12_crossref_pass.txt all
```
The resulting files will be placed in your output folder - all.events and all.quant.

#### VALIDATOR PARSE
##### Requirements
Spreadsheet-ParseExcel at least version 0.59. This module has several pre-requisites of its own. These files are located in the folder code/perl_modules_required.
##### Usage
* Create a file called ValidatorList in the supporting_files folder contains:
 * First Line: Path to excel output file from Validator. This path should be relative to the data folder. 
 * Second line, list he name of the excel output file.

```
perl ValidatorParse.pl ../supporting_data/ValidatorList
```
Three files will be created:
* AllMatch
* HLRecover
* MascotMatch

### Usage - Toolchain

The software stack can also be run with all of the file parameters provided on the command line. To see the command line options run:

```python validator_cli.py --help```

By default validator_cli.py only Validator runs. If you include a CPAS/Peptide Prophet file as a parameter validator_cli.py will also run Crossref, Quant, and Validator Parse



 