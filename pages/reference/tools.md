# Bioinformatics Tools

Summary of purpose & usage of bioinformatics tools seen in class.

## IDBA-UD
**Assembles reads into contigs**

Example usage (replace names in brackets with your own names and remove brackets):

```
idba_ud -r [reads file] -o [output folder name]
```

* The `-r` gives it the "path" (directions) to the reads for your reads. `../` means it is in the directory outside of the one you're in.
* The `-o` flag tells the program what you want the output directory to be called.


## Prokka
**Finds ORFs in assembled contigs and annotates them**

Example usage (replace names in brackets with your own names and remove brackets):
```
prokka [assembled contigs fasta file] --outdir [name of output directory]
```

To get a summary of how many ORFs were assigned to different COG categories based on your Prokka data, do this:
```
get_ORF_COG_categories.py [name of your Prokka tsv file]
```

## BLAST
**Find a match to a sequence in your own sequence file, not the NCBI database**

Example usage (replace names in brackets with your own names and remove brackets):

First, make a BLAST database:
```
makeblastdb -in [example sequence file] -dbtype [prot or nucl]
```
* -in gives the name of the file you want to BLAST against.
* for dbtype, use `prot` for an amino acid file and `nucl` for a nucleotide file.

Then, do your BLAST.

### blastp
```
blastp -query [example query file] -db [example database file] -evalue 1e-05 -outfmt 6 -out [example_output.blastp]
```
- `blastp` does a protein vs protein blast. (other choices are `blastn`, `blastx`, `tblastn`, `tblastx`.)
- `-query` defines your blast query-- in this case, the Pfam seed sequences for the CRISPR RAMP proteins.
- `-db` defines your database-- in this case, the toy assembly ORFs.
- `-evalue` defines your maximum e-value, in this case 1x10-5
- `-outfmt` defines the output format of your blast results. option 6 is common; you can check out https://www.ncbi.nlm.nih.gov/books/NBK279675/ for other options.
- `-out` defines the name of your output file. I like to title mine with the general format `query_vs_db.blastp` or something similar.

### tblastn

```
tblastn -query [example query file] -db [example database file] -evalue 1e-05 -outfmt 6 -out [example_output.tblastn]
```
- `tblastn` does a protein vs translated nucleotide blast. (other choices are `blastn`, `blastx`, `blastp`, `tblastx`.)
- `-query` defines your blast query-- in this case, the Pfam seed sequences for the CRISPR RAMP proteins.
- `-db` defines your database-- in this case, the toy assembly ORFs.
- `-evalue` defines your maximum e-value, in this case 1x10-5
- `-outfmt` defines the output format of your blast results. option 6 is common; you can check out https://www.ncbi.nlm.nih.gov/books/NBK279675/ for other options.
- `-out` defines the name of your output file. I like to title mine with the general format `query_vs_db.blastp` or something similar.


## Making a tree
**First make an alignment with muscle, then make a tree with RAxML**

First, use nano or some other text editing tool to make a fasta file with your sequences of interest. The following is an example of how to make a tree using amino acid (not nucleotide) sequences.

Example usage (replace names in brackets with your own names and remove brackets):

```
muscle -in [example_sequence_file.fasta] --out [example_sequence_file_aligned.afa]
convert_afa_to_phy.py [example_sequence_file_aligned.afa]
raxmlHPC-PTHREADS-AVX -f a -# 20 -m PROTGAMMAAUTO -p 12345 -x 12345 -s [example_sequence_file_aligned.phy] -n [example_tree_name.tree] -T 4
```
Open the RAxML_bipartitions.example_tree_name.tree file in ITOL.


## Mapping
**Map short reads against a reference**

Example usage (replace names in brackets with your own names and remove brackets):

```
bowtie2-build [reference fasta file] [reference_file.btindex]
bowtie2 -x [reference_file.btindex] -f -U [reads_to_map.fasta] -S [output_file.sam]
samtools view -bS [output_file.sam] > [output_file.bam]
samtools sort [output_file.bam] -o [output_file_sorted.bam]
samtools index [output_file_sorted.bam]
```

To get ORF coverage:

```
make_bed_file_from_gff_file_prokka.py [your gff file from prokka]
samtools bedcov [your bed file] [your sorted bam file] > [ an output file that ends in _ORF_coverage.txt]
```

## anvi'o
**Make bins from metagenomic assemblies and BAM files**

The protocol for week 6 laid out anvi'o step by step. The [anvi'o tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) is another good source of information.


## mothur
**Analyze 16S data to find taxonomy, OTUS**

The protocol for week 7 laid out mothur step by step. The [mothur wiki](https://www.mothur.org/wiki/MiSeq_SOP) is another good source of information.


# Legacy Items

## Prodigal

**Finds ORFs in assembly file**

Example usage:

```bash
mkdir ORF_finding
cd ORF_finding
prodigal -i ../toy_assembly/toy_dataset_assembly_subsample.fa -o toy_assembly_ORFs.gbk -a toy_assembly_ORFs.faa -p single
```

* The `-i` flag gives the input file, which is the assembly you just made.
* The `-o` flag gives the output file in Genbank format
* The ‘-a” flag gives the output file in fasta format
* The `-p` flag states which procedure you’re using: whether this is a single genome or a metagenome. This toy dataset is > a single genome so we are using –p single, but for your project dataset, you will use –p meta.


## Interproscan

**ORF Fasta File -> Annotated TSV**

Used to efficiently and effectively annotate proteins. Compares your open reading frames against several protein databases and looks for protein "signatures," or regions that are highly conserved among proteins, and uses that to annotate your open reading frames. It will do this for every single open reading frame in your dataset, if it can find a match.

Example usage:

```bash
interproscan.sh -i toy_assembly_ORFs.noasterisks.faa -f tsv
```

* The `-i` flag gives the input file, which is your file with ORF sequences identified by Prodigal, with the asterisks removed.
* The `-f` flag tells Interproscan that the format you want is a tab-separated vars, or “tsv,” file.
