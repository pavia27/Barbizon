### Extracts genomic sequence composition and build a genetic landscape of genetic information

By providing a folder containing either metagenome assembled genomes (MAGs) or the contigs generated from an assembly this program will output a file in TSV format summarizing GC%, tetranucelotide frequency, codon frequency, depth (optional), and the taxanomic group (optional) matched to each contig.

***

### Easy installation:

```
git clone https://github.com/pavia27/Barbizon.git
cd Barbizon
bash setup.sh
source activate barbizon
```

*Do not worry about the dependencies after conda installation. Just enter source deactivate when finished using the program.

***

###  Usage

Print out the help screen :

```
./barbizon.py -h
```

For quick usage to generate a tsv of GC%, tetranucleotide identity, and codon usage :

```
./barbizon.py -b MAGS/ -o test_MAG_set
```
To include depth information :

```
./barbizon.py -b MAGS/ -o test_MAG_set --bam optional.bam.sorted
./barbizon.py -b MAGS/ -o test_MAG_set --bam optional_1.bam.sorted,optional_2.bam.sorted,optional_3.bam.sorted
```

To include taxonomic information :
Please provide a reference protein dataset for use by CAT. Below are commands for downloading the database:
```
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz
tar -xvzf CAT_prepare_20210107.tar.gz
```
Move CAT_prepare_20210107 directory to the Barbizon directory 

```
./barbizon.py -b example_MAGS/ -o test_MAG_set --bam optional.bam.sorted -t 10 -c TRUE
```
