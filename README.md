## acCRISPR: an activity-correction method for improving CRISPR screen accuracy

### About 

### Installation
1. `git clone https://github.com/ianwheeldon/acCRISPR.git`
2. `cd acCRISPR`
3. `pip install .`

### Usage

### Example run

Navigate to directory containing source code files: `cd src/acCRISPR`

Copy test input files into this directory: `cp ../../example_data/*.tab .`

Finally, run acCRISPR using the following command:
```bash
./run_acCRISPR.py --counts pH3_counts_final.tab --replicate_info pH3_rep_file.tab --cov 6 --significance 2-tailed --output_prefix test_run
```
