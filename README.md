# BRCA1_ecDNA_custom_scripts
Custom scripts used in BRCA1_ecDNA project

## script for filtering cnvkit results to generate AA seeds
trace_ampliconsuite_cnv_to_seeds.py

This script filters CNVkit results to generate AmpliconArchitect (AA) seeds and reports the specific filtering reason for each discarded CNV record. It provides the flexibility to bypass the standard arm-median-based threshold and BAM coverage-ratio check originally implemented in the AmpliconSuite-pipeline.

### Requirements:
AmpliconSuite-pipeline (https://github.com/AmpliconSuite/AmpliconSuite-pipeline) version 1.5.0

Python 3.10.18

### Installation:
No installation needed.

### Example Usage:
Run with bam coverage check:
```
python3 trace_ampliconsuite_cnv_to_seeds.py --cns BZ2213695K_WT.cs.rmdup.cns --ref hg19 --aa_data_repo /path/to/AA_data_repo --gain 4.0 --cnsize_min 50000 -s BZ2213695K_WT --outdir ./BZ2213695K_WT BZ2213695K_WT --bam BZ2213695K_WT.cs.rmdup.bam 1> stdout.txt 2> stderr.log
```

Run without bam coverage check and disable extra chromosome-arm-cnv based threshold:
```
python3 trace_ampliconsuite_cnv_to_seeds.py --cns BZ2213695K_WT.cs.rmdup.cns --ref hg19 --aa_data_repo /path/to/AA_data_repo --gain 4.0 --cnsize_min 50000 -s BZ2213695K_WT --outdir ./BZ2213695K_WT --disable_stage1_med_ccg_thresh --disable_stage2_bam_coverage_check 1> stdout.txt 2> stderr.log
```

The .cns and .bam files are generated in AmpliconSuite pipeline by default.

By default, this scirpts using --bam to do coverage check as the same as AmpliconSuite. This could be skipped by using --disable_stage2_bam_coverage_check. It typically takes few seconds to one minite to run without --bam, and it may take few minutes when using --bam, depending on the size of the bam file.

The results containing _AA_CNV_SEEDS.bed file which could be used as the input of --cnv_bed in AmpliconSuite-pipeline, and 2 other files .trace_SEED_stage1.tsv / .trace_SEED_stage2.tsv containing reasons of filtering.

