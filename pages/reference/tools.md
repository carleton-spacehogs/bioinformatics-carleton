# Bioinformatics Tools

Summary of purpose & usage of bioinformatics tools seen in class.

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
