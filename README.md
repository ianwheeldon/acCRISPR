# acCRISPR: an activity-correction method for improving CRISPR screen accuracy

acCRISPR is a tool developed in Python 3 for analyzing genome-wide CRISPR screens to identify essential genes for growth as well as gain- and loss-of-function hits for stress tolerance by removing low-activity sgRNA based on an activity cutoff. acCRISPR optimizes library activity and coverage to improve hit calling and screen accuracy. 

## System requirements
acCRISPR can be implemented on any desktop computer running macOS or Linux. The tool has been tested on Linux: Ubuntu 20.04.

## Dependencies
acCRISPR requires Python version 3.6 or higher and the following Python libraries to be installed:
- `numpy (>=1.18.1)`
- `scipy (>=1.4.1)`
- `matplotlib (>=3.4.2)`
- `statsmodels (>=0.12.0)`

## Installation
Steps:
1. Clone this repository: `git clone https://github.com/ianwheeldon/acCRISPR.git`
2. Navigate to acCRISPR directory: `cd acCRISPR`
3. Install the Python package: `pip install .`

Provided that all dependencies have already been installed, acCRISPR installation usually takes only about 1-2 mins. After executing the above steps, installation of the package can be confirmed by launching Python3 on the terminal and running `import acCRISPR`. If this results in no error, the package was installed correctly.

## Running acCRISPR on an example dataset

*(Note: For the purposes of acCRISPR implementation, 'FS' denotes fitness score when analyzing growth screen data and tolerance score when analyzing stress tolerance screen data)*

Navigate to the directory containing source code: `cd src/acCRISPR`

Copy test input files into this directory: `cp ../../example_data/*.tab .`

The example dataset used here is the pH 3 dataset from tolerance screening experiments in *Yarrowia lipolytica* using CRISPR-Cas9. acCRISPR can be run for the original/uncorrected sgRNA library using the following command (run time ~ 7-10 mins):
```
python3 run_acCRISPR.py --counts pH3_counts_final.tab --replicate_info pH3_rep_file.tab --cov 6 --significance 2-tailed --output_prefix pH3_no_cutoff
```
The input files used in the above example can be found in the `example_data` directory.

- `--counts` is used to specify the file containing raw sgRNA counts for all replicates of each sample used in the screening experiment. The first and second columns specify the guide identifiers and associated gene names respectively, and the subsequent columns contain information on raw counts for each sgRNA in each sample replicate.
- `--replicate_info` contains information for mapping replicates to samples. The names of replicates in this file and their order should exactly match the column headers (for replicates) and their order (from left to right) in the file specified by `--counts`. If different controls are used for determining CS & FS, this file should contain 4 unique sample names - `Control_CS`, `Treatment_CS`, `Control_FS` & `Treatment_FS`. However, if a common control sample is used for both CS & FS estimation, only 3 unique sample names are needed - `Control`, `Treatment_CS` & `Treatment_FS`.
- `--cov` is used to provide the coverage (avg. no. of sgRNA per gene) of the original library, rounded off to the nearest integer. In the given example, the original library contains 46395 unique sgRNA targeting 7854 genes in the genome. Hence, the library coverage would be 5.91, which can be rounded off to 6.

For brief information on all input parameters to acCRISPR: `python3 run_acCRISPR.py -h`

If CS and FS of sgRNA have been pre-calculated, CS & FS (i.e., log2-fc) calculation can be skipped and the file containing CS & FS values can be provided directly as input to acCRISPR. In order to skip log2-fc calculation, the parameter `--skip_log2fc_calc` needs to be `True`.
```
python3 run_acCRISPR.py --skip_log2fc_calc --CS_FS_file CS_FS_values.tab --cov 6 --cutoff 5.0 --significance 2-tailed --output_prefix pH3_5.0
```
In the above command, the filename specified for the parameter `--CS_FS_file` (i.e., `CS_FS_values.tab`) should contain CS and FS values of sgRNA. This file should be in the same format as `pH3_no_cutoff_guide_CS_FS.tab` found in the directory `example_no_cutoff`.

To run acCRISPR on the example pH 3 dataset with a CS-corrected library at a threshold of 5.0 & providing raw sgRNA counts as inputl, the `--cutoff` parameter should be set to 5.0, as in the command below:
```
python3 run_acCRISPR.py --counts pH3_counts_final.tab --replicate_info pH3_rep_file.tab --cov 6 --cutoff 5.0 --significance 2-tailed --output_prefix pH3_5.0
```
Alternatively, acCRISPR can be implemented for a corrected library by skipping CS & FS calculation from raw counts and providing CS & FS of sgRNA from the original/uncorrected library as input instead of the raw counts file. For the given example dataset, once acCRISPR has been implemented with the uncorrected library, it generates a file named `pH3_no_cutoff_guide_CS_FS.tab`, which can be provided as input to generate results at a threshold of 5.0 using the `--CS_FS_file` parameter. In addition, the parameter `--skip_log2fc_calc` needs to be `True`.
```
python3 run_acCRISPR.py --skip_log2fc_calc --CS_FS_file pH3_no_cutoff_guide_CS_FS.tab --cov 6 --cutoff 5.0 --significance 2-tailed --output_prefix pH3_5.0
```

acCRISPR provides the following output files:

- A TAB delimited file containing CS and FS of sgRNA in the corrected library (if a CS threshold is provided & `--skip_log2fc_calc` is `True`) or uncorrected library (if no CS threshold is provided). The file contains 4 columns - sgRNA ID, associated gene name, CS & FS. (If a CS threshold is provided but `--skip_log2fc_calc` is `False`, 2 separate files containing CS & FS of uncorrected and corrected libraries are output.)
- A TAB delimited file containing (uncorrected or corrected, depending on whether a threshold is provided) FS of all genes covered by the library. The file contains 2 columns - gene name & gene FS.
- A TAB delimited file containing raw and corrected p-values of all genes. The file contains 4 columns - gene name, gene FS, raw p-value and corrected p-value.*
- A PNG file showing the null distribution plot and corresponding null distribution parameters (mean & S.D.), as well as the number of sgRNA used per pseudogene to create this null distribution.*
- A text file containing information about the number of significant genes, average coverage of the uncorrected/corrected library, and the value of ac-coefficient.*

<b><i>*Since the algorithm for hit identification by acCRISPR involves random sampling, the null distribution S.D. and hence, the p-values of genes and no. of significant genes would vary slightly with each run of the tool.</i></b>

Output files obtained from runs on the example dataset with the uncorrected and corrected (threshold = 5.0) libraries can be found in `example_no_cutoff` and `example_cutoff_5.0` directories respectively.
