## acCRISPR: an activity-correction method for improving CRISPR screen accuracy

### About 

### Installation
1. `git clone `
2. `cd acCRISPR`
3. `pip install .`

### Usage

### Example run
1. `cd src/acCRISPR`
2. `cp ../../example_data/*.tab .`
3. `./run_acCRISPR.py --counts pH3_counts_final.tab --replicate_info pH3_rep_file.tab --cov 6 --significance 2-tailed --output_prefix test_run`
