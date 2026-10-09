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
- `ROW_NUMBERS` - number the rows of the plot (the numbers match the index of the `lrRNAseq_plot_df` table, whose `TrNames` column lists the transcripts in each row).
- `PLOT_PATH` - file to save the plot to (`.png`, `.pdf` or `.svg`); leave empty to only display it.

# Reading the plot
- **Red blocks:** regions of the genome where a part of a transcript aligned.
- **Coloured lines:** the extent of each transcript; each transcript in a row has its own colour.
- **Gold outline / gold fill:** the aligned part overlaps one of the transcript's ORFs / its largest ORF.
- **Numbers:** the position of each aligned region within its transcript (0 = nearest the 5' end). Along the genome, `0 1 2 3 ...` is a gene on the forward strand, `... 3 2 1 0` a gene on the reverse strand, and out-of-order numbers show rearrangements or repeats inside the transcript. Regions where the same part of a transcript aligned more than once share a number.
- **Background:** colinear blocks of regions - pink where the numbers increase along the genome (or a single region), blue where they decrease. A transcript split into several blocks has a rearrangement or repeat.

The numbering assumes transcripts are oriented 5'→3', as full-length (FLNC/CCS) reads normally are.

[Full documentation](Docs/overview.md)


<br />
