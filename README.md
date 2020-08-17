# human-atac
## Cell barcode and demultiplexing for sci-ATAC-seq3 runs

### Overview
We recommend first demultiplexing using BCL2FASTQ with a placeholder samplesheet and then running a second script to split reads by sample and assign barcodes. From there, generic single-cell ATAC-seq pipelines like [SnapATAC](https://github.com/r3fang/SnapATAC) should be able to process the data further into count matrices, etc.

#### BCL2FASTQ
```
bcl2fastq --runfolder-dir <PATH_TO_BCL_DIRECTORY> \
               -o <OUTPUT_DIRECTORY> \
               --interop-dir <OUTPUT_DIRECTORY> \
               -r 6 -w 6 \
               --sample-sheet <PLACEHOLDER_FILE> \
               --ignore-missing-positions \
               --ignore-missing-controls \
               --ignore-missing-filter \
               --ignore-missing-bcls
```

Note that in this case we provide a dummy samplesheet (`samplesheet_fake.txt`) for the --sample-sheet argument to bcl2fastq that looks like:
```
[DATA],,
Sample_ID,Sample_Name,index,index2
fake,fake,NNNNNNNNNNNNNNNNNNNN,NNNNNNNNNNNNNNNNNNNN
```

Note that the specified indices must be 20bp in length (and therefore have N(20) in file above) for 3-level sci-ATAC-seq experiments when following our protocol.

This will split the fastq files out by lane, which can be useful for parallelization of downstream steps, but bcl2fastq has a `--no-lane-splitting` option that will output to a single file if you prefer.

#### Barcode assignment and demultiplexing
We provide `barcode_correct_sciatac.py` for demultiplexing sci-ATAC-seq3 data by sample and assigning cell barcodes.

##### Sample sheet
For the following step you must provide a sample sheet that specifies the combinations of ligation and PCR indices that define each sample. The sample sheet should specify the sample_id for each sample in your experiment and the ranges of the barcodes used for the corresponding well positions. The basic format is: [N7 ligation barcode]:[P7 PCR barcode]:[P5 PCR barcode]:[N5 ligation barcode]. Note that the N7 and P7 barcodes are numbered row-wise (well A1 = 1, well A2 = 2, …, well B1 = 13 etc.) whereas N5 and P5 barcodes are numbered column-wise (well A1 = 1, well B1 = 2, …, well A2 = 9 etc.). You can specify the barcodes as ranges or list them individually separated by commas.

```
sample_id ranges
sample_id_1	N7:P7:P5:N5 
sample_id_2	N7:P7:P5:N5 
```

Sample sheets for this study are provided as `samplesheet_batch1_atlas.txt`, `samplesheet_batch2_atlas.txt` and `samplesheet_batch3_atlas.txt`. 

In this case, 24 samples were laid out across one 96-well plate as described in [protocols.io](https://dx.doi.org/10.17504/protocols.io.be8mjhu6). 

Since samples are pooled after the first round of indexing (N5 ligation), only the N5 ligation indices will differ across samples. For the remaining sets of indices include all used barcodes for every sample. See `samplesheets` and below for examples.

e.g.

```
sample_id ranges
sample_1	1-384:1-48:1-32:1,9,17,25,97,105,113,121,193,201,209,217,289,297,305,313 
sample_2	1-384:1-48:1-32:2,10,18,26,98,106,114,122,194,202,210,218,290,298,306,314
```
Here, sample_1 was layed out in wells A1-A4 for the initial round of indexing, and sample 2 was layed out in wells B1-B4, as specified by the column-ordered N5 indices.


##### Script
```
python barcode_correct_sciatac.py \
  --samplesheet <SAMPLESHEET_PATH> \
  -1 <(zcat <READ1_FASTQ_GZIPPED>) \
  -2 <(zcat <READ2_FASTQ_GZIPPED>) \
  --output1_prefix r1. \
  --output2_prefix r2. \
  --stats_out <OUTPUT_STATS_JSON> \
  --wells_384
```

Here we provide the samplesheet described above (for this script and not the placeholder for bcl2fastq), the two fastq files for R1 and R2 via process substitution (`<()`), a prefix for the output `r1` and `r2` files (defaults are the values specified above), a JSON file that will be written to provide output stats about barcode correction, a flag to specify that the barcode set we are using is the 384 well set used for fetal samples in the paper. You may optionally provide the flag `-X` if the run was on a NextSeq or similar where some of the barcodes will need to be reverse complemented and retrieved differently.

The script could be modified to accept file names directly rather than the contents via process substitution -- we found this implementation to be somewhat faster and more flexible for executing small tests (e.g. could modify to use `-1 <(zcat <READ1_FASTQ_GZIPPED> | head -10000)` to run on small subset).

In practice we also ran this script using [pypy](https://www.pypy.org/) internally and found it to be much faster.

See `python barcode_correct_sciatac.py --help` for more script documentation.

##### Outputs
We output a JSON file that looks roughly like the following:
```
{
    "pcr_i5": 0.96239,
    "pcr_i7": 0.93147,
    "tagmentation_i5": 0.9806,
    "tagmentation_i7": 0.95688,
    "all_barcodes": 0.89864,
    "total_input_reads": 100000,
    "total_not_specified_in_samplesheet": 0
}
```

Where we note the total number of input reads, the fraction with valid barcode assignments for each portion of the combinatorial cell barcode, and the fraction with a valid overall cell barcode assignment (valid assignment for all parts; `all_barcodes`).

It will also output a pair of fastq files (r1 and r2) for each sample noted in the provided samplesheet. For example, if your samplesheet includes a sample named `sample_70_lung`, it would output `r1.sample_70_lung.fastq.gz` and `r2.sample_70_lung.fastq.gz` given the command above. These files are output uncompressed and then gzipped at the end, which we found to be faster than directly writing gzipped output at the expense of more disk usage. You could modify the script to directly write gzipped output, use alternate compression tools like `pigz`, or parallelize the compression step of this script as needed.

##### Script testing
In addition, we provide small test fastqs (`Undetermined_R1_for_test.fastq.gz` and `Undetermined_R2_for_test.fastq.gz`) which allow you to test the demux script on a smaller subset of 100000 reads (the appropriate samplesheet for these files is provided as `samplesheet_test_demux.txt`).

#### Downstream analysis
As mentioned above, there are already well-supported tools for processing scATAC-seq data such as [SnapATAC](https://github.com/r3fang/SnapATAC). After demultiplexing we recommend that you utilize an available tool to generate count matrices, etc.

Tools like SnapATAC as well as [Signac](https://satijalab.org/signac/) (which integrates very nicely with existing Seurat workflows) also now offer excellent support for downstream analysis like visualization, dimensionality reduction, clustering, annotation, etc.

