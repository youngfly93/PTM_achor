# Results Directory

This directory contains the output files from the PTM-HLA anchor coupling analysis pipeline.

## Expected Output Files

After running the pipeline, you will find the following files here:

### Core Analysis Results
- `tumor_anchor_mod.tsv` - Tumor sample modification-anchor coupling records
- `normal_anchor_mod.tsv` - Normal sample modification-anchor coupling records  
- `all_anchor_mod.tsv` - Combined modification-anchor coupling records

### Statistical Analysis
- `anchor_mod_enrichment.tsv` - Modification enrichment analysis at anchor positions
- `anchor_mod_group_diff.csv` - Tumor vs Normal Fisher exact test results

### Reports
- `tumor_vs_normal_report.txt` - Detailed comparison analysis report
- `pipeline_summary.txt` - Pipeline execution summary
- `Phospho_anchor_violin.txt` - Phosphorylation modification visualization (if data available)

## File Formats

### TSV Files
Tab-separated values with headers containing:
- Sequence information
- HLA allele predictions  
- Modification details
- Anchor position annotations
- Statistical scores

### CSV Files
Comma-separated values with statistical test results and effect sizes.

### TXT Files
Human-readable reports with analysis summaries and interpretations.

## Usage

These files are generated automatically when running:

```bash
python3 run_pipeline.py
```

For test runs, use:
```bash
python3 run_pipeline.py --test
```