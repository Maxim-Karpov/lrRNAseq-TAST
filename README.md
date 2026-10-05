# Overview
LrRNAseq TAST processes BLAST alignments of full length long read RNAseq transcripts to a relatively short genomic region (e.g. 100kb) and facilitates recognition of gene splicing.

![Screenshot](318A10_cov25_wORFs_wSynteny_wSpacing_wSplicing_normal_intron_filter_highlights2.png)

The **lrRNAseq_TAST_visualiser.ipynb** is the main notebook used for visualisation of collective BLAST alignments (total alignment). The notebook contains detailed instructions of use within.

The **bioinformatic_functions.ipynb** is required to be imported in the main file to call multiple functions from.

The rest of the files serve as example data to try out.

# Requirements
- Python 3.9 or newer, with Jupyter Notebook or Jupyter Lab
- Python libraries: Pandas, Matplotlib, Numpy, Biopython. Install them with `pip install -r requirements.txt`
- BLAST, outputting in -outfmt '6 qseqid sseqid slen qlen qstart qend sstart send length mismatch gapopen pident evalue bitscore' format (no header line).

# Main settings
All settings are written in capital letters in the visualiser notebook. The most important ones are:
- `STRAND_MODE` - which alignment strand(s) to keep for each transcript. `"dominant"` (default) keeps only the strand with the highest total bitscore, so sense and antisense hits are not merged into one gene model. `"both"` reproduces the behaviour of version 1.0.
- `MAX_INTRON_LENGTH` and `MIN_COVERAGE` - transcript filters.
- `BARCODE_LEN` - number of barcode bases trimmed from both ends of each transcript before ORF finding. ORF positions are automatically shifted back so they match the BLAST coordinates.
- `MIN_GAP` - minimal genomic distance between transcripts that share a plot row.

[Full documentation](Docs/overview.md)


<br />
