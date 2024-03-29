# Week 5: Mapping with bowtie2

#### 1. Log in to the remote server
Boot your computer as a Mac and use the Terminal to ssh in to baross.

## Mapping with a toy dataset

#### 2. Make new directory
We're going to start by mapping the sequencing reads from a genome sequence of a type of archaeon (*Sulfolobus acidocaldarius*) against a scaffold from a very closely related species.

The sequencing reads from from one strain of *Sulfolobus acidocaldarius*, and the reference sequence that they are mapping to is from a very closely related strain of *Sulfolobus acidocaldarius*.

In your toy dataset directory, make a new directory called “mapping,” then change into that directory.
```
cd ~/toy_dataset_directory
mkdir mapping
cd mapping
```
#### 3. Copy data
Copy this week's toy datasets into your directory. You're going to copy:
- the reference sequence (toy_dataset_contig_for_mapping.fasta)
- the reads that you will map to the reference (toy_dataset_reads_for_mapping.fasta)
```
cp /usr/local/data/toy_datasets/toy_dataset_reads_for_mapping.fasta .
cp /usr/local/data/toy_datasets/toy_dataset_contig_for_mapping.fasta .
```

#### 4. Build index
You're going to be mapping short reads to a longer reference sequence. The first thing you have to do is prepare an **index** of your reference so that the mapping software can map to it.

```
bowtie2-build toy_dataset_contig_for_mapping.fasta toy_dataset_contig_for_mapping.btindex
```
- `bowtie2-build` is the program that indexes your reference.
- The first argument gives the reference dataset name.
- The second argument provides the name you want to give to the index.


#### 5. Map!
Now, map! This will take a few steps.

First, you make what is called a SAM file. It's a human-readable version of a BAM file, which we read about in the Zimmer "Game of Genomes" articles.
- `bowtie2` is the name of the mapping program.
- `-x` is the flag that provides the name of the index you just made.
- `-f` means that the reads you are mapping are in fasta, not fastq, format.
- `-U` means that the reads are not paired. (They aren't in this dataset.)
- `-S` provides the name of your output file, which is in SAM format.

```
bowtie2 -x toy_dataset_contig_for_mapping.btindex -f -U toy_dataset_reads_for_mapping.fasta -S toy_dataset_mapped_species1.sam
```

#### 6. Look at SAM file output
If you look at the output file with `less`, you can see that it is human-readable (sort of). It tells you exactly which reads mapped, and where they mapped on the reference, and what the differences were between the reference and the mapped reads. This can sometimes be useful if you want to parse it yourself with your own scripts-- but there's a whole suite of tools in a package called `samtools` that we'll rely on to do that next.

#### 7. samtools
Now you will use a package called `samtools` to convert the SAM file into a non-human-readable BAM file. Carl Zimmer spent a lot of time trying to obtain his BAM files in that "Game of Genomes" article-- now you get to make one yourself! (But, uh, not of your own genome, obviously...)

```
samtools view -bS toy_dataset_mapped_species1.sam > toy_dataset_mapped_species1.bam
```

- ` samtools` is a package used to manipulate and work with mapping files. samtools view is one program within the whole samtools package.
- The flag `-bS` is not BS! It tells samtools to convert a sam file to a bam file. (Bioinformatics jokes = still not very funny.)


#### 8. samtools sort
And because this is so fun, we get to do some more bookkeeping. Sort your bam file so that later programs have an easier time parsing it:

```
samtools sort toy_dataset_mapped_species1.bam -o toy_dataset_mapped_species1_sorted.bam
```

- `samtools sort` is the name of the program used for sorting
- The first argument provides the name of the bam file you want to sort
- The `-o` flag gives the name of the output file you want.

#### 9. Index the reference with samtools
In order to visualize your mapping, you have to **index your reference**. Because indexing is SO MUCH FUN. This time we are indexing with samtools instead of bowtie2.

(**NOTE!!**: you only need to do this step if you're going to visualize your mapping with IGV, as we're about to do now. In the future, if you don't intend to visualize your mapping, then you don't need to bother with this step.)

```
samtools faidx toy_dataset_contig_for_mapping.fasta
```

- `samtools faidx` is the name of the program that indexes the reference.
- The first argument provides the name of the index, which should be your reference file.

#### 10. Index the bam file
Almost there! Now you index the bam file that you just made because WE LOVE INDEXING.

```
samtools index toy_dataset_mapped_species1_sorted.bam
```

- `samtools index` is the name of the program that indexes the bam files.
- The first argument provides the name of a sorted bam file.

#### 11. Copy to local computer
Now we're going to visualize this. Copy the entire "mapping" folder over to your local computer using scp (or Filezilla if you prefer.)

```
scp -r username@baross.its.carleton.edu:~/mapping/ .
```
*Remember, to use scp, you should open a new Terminal window that is NOT logged in to baross. The above command copies the folder at toy_dataset_directory/mapping to your local computer in whatever folder you happen to be in. You could tell the computer to put it on the Desktop or wherever you like by simply providing the path instead of a period in the second part of the command.*

#### 12. Visualize in IGV
Find the IGV Viewer on your local computer and open it. **First,** click 'Genomes' --> 'Load Genome from File' and find your reference file (`toy_dataset_contig_for_mapping.fasta`). **Next,** click 'File' --> 'Load from File' and open your sorted bam file (`toy_dataset_mapped_species1_sorted.bam`). You should be able to visualize the mapping. Along the top, you'll see the coordinates of your reference sequence. Below that, you'll see a graph showing the coverage of each base pair along your reference sequence. Below that, you'll see each read mapped to each position. The arrows indicate the direction of the read; white reads are reads that mapped to two different locations in your reference. Single nucleotide variants in the reads are marked with colored letters; insertions are marked with a purple bracket, and deletions are marked with a horizontal black line. More information can be found at the link below.

Link to [IGV viewer](http://software.broadinstitute.org/software/igv/AlignmentData)

#### 13. Compare mappings
We're going to compare and contrast this mapping with another one. Now we're use the sequencing reads from a third very closely related strain of *Sulfolobus acidocaldarius*, and we're going to map those reads to the original reference sequence so that we can compare the mapping.

All of the commands are listed below. First, you will copy the second file to your directory (see below). You have already indexed the reference file. Then you map, convert SAM to BAM, sort, and then index.

```
cd ~/mapping/
cp /usr/local/data/toy_datasets/toy_dataset_reads_for_mapping_species2.fasta .
bowtie2 -x toy_dataset_contig_for_mapping.btindex -f -U toy_dataset_reads_for_mapping_species2.fasta -S toy_dataset_mapped_species2.sam
samtools view -bS toy_dataset_mapped_species2.sam > toy_dataset_mapped_species2.bam
samtools sort toy_dataset_mapped_species2.bam -o toy_dataset_mapped_species2_sorted.bam
samtools index toy_dataset_mapped_species2_sorted.bam
```

#### 14. Visualize new mapping
Copy the new data files to your local computer using FileZilla or scp like in step 11. Visualize both of them in IGV viewer. Since you have already loaded the reference file and reads from your first mapping, all you have to do is click 'File' --> 'Load from File' and click on `toy_dataset_mapped_species2_sorted.bam`. You should be able to see them side by side.

**Pause here. Check in with your group. Together, discuss and answer these questions on a shared Google Doc that you share with Rika.**

**Check for understanding:**

Q1. Describe the large-scale differences between the mapped reads from species 1 and species 2, and explain what this mapping tells us about the relative genome structure of the two genomes that we mapped. If we compared this genomic region in a dot plot, what would it look like?

Q2. Do you see evidence of misassemblies or deletions in the reference? What does that evidence look like?



## Mapping your project datasets
Now we're going to map your project datasets. Remember that these are metagenomes, not a genome, so the data will be a bit more complex to interpret.

We're going to map your *reads* against your *assembled contigs*. Why would we do this, you ask? A few reasons:
- to look for single nucleotide variants in specific genes.
- to quantify the relative abundances of different genes, and determine whether specific genes have better coverage than others.
- to quantify the relative abundances of specific taxa, and determine whether specific taxa are more abundant than others.

**As you consider this, discuss the following question with your group:**

**Check for understanding:**

 Q3. If you wanted to quantify the relative abundances of specific genes in your sample, why couldn't you simply count the number of times your gene appears in your assembly?

#### 15. Map to your project datasets
Change directory into your project dataset directory folder.

```
cd ~/project_directory
```

We're going to map your raw reads against your assembled contigs (not your ORFs). Make sure you know where your project assembly is and where your raw reads are. Follow the instructions to map your raw reads back to your assembled contigs. For example, if you were mapping the reads in the dataset `ERR599166_1mill_sample.fasta` against an assembly called `ERR599166_assembly_reformatted.fa`, you might do something like this (below). 

One more thing. Your project datasets are very large-- you each have 10 million reads. So mapping will take longer than for the toy datasets. So we're going to add an extra flag (-p) to the mapping step to tell the computer to use more than one CPU so that this goes faster. You'll each use 4 CPUs for this process. The mapping may still take about 5-10 minutes, so have patience!

An example set of commands is shown below. **Remember to replace the dataset names here with your own project datasets!**

```
bowtie2-build ERR599166_assembly_reformatted.fa ERR599166_assembly_reformatted.btindex
bowtie2 -x ERR599166_assembly_reformatted.btindex -f -U ERR599166_sample.fasta -S ERR599166_mapped.sam -p 4
samtools view -bS ERR599166_mapped.sam > ERR599166_mapped.bam
samtools sort ERR599166_mapped.bam -o ERR599166_mapped_sorted.bam
samtools faidx ERR599166_assembly_reformatted.fa
samtools index ERR599166_mapped_sorted.bam
```

*Note that this command assumes that your raw reads and your assembly are in the same directory you're in. If they are not, you will need to either copy them over or use the correct path in your commands. (A reminder: the path simply gives directions to the computer for where to find a file.) For example, if you are in a mapping directory, your reads file is one directory up in the file hierarchy, and your assembled reads are in your assembly file, you might have to type something like this:*

`bowtie2 -x ../assembly/ERR599166_assembled.btindex -f -U ../ERR599166_1mill_sample.fasta -S ERR599166_mapped.sam`

#### 16. Visualize
When you visualize this in IGV, remember that you have multiple contigs. So you have to click the drop-down menu at the top and choose which contig you wish to visualize.

#### 17. Check for understanding

**With your table, discuss and answer the following and put the answer in the shared Google Doc:**

 Q4. Do you see evidence of single nucleotide variants? Biologically speaking, what does this indicate? (Keep in mind that you have mapped metagenomic reads from a whole microbial community against a consensus assembly-- this is not reads from an individual vs an individual's reference assembly.)

#### 17. Calculating coverage- generate bed file
You were able to visualize the mappings in IGV, but sometimes you just want to have a number: for example, you might want to know the average coverage across a specific gene, and compare that to the average coverage of another gene in order to compare their relative abundances in the sample. So, next we're going to calculate gene coverages based on your mapping.

First we're going to make what's called a bed file. We will use it to find the average coverage of every single open reading frame in your dataset. Please make sure that your contigs have names that are something like `c_00000000001` and your ORFs have names that are something like `c_00000000001_1`.

```
make_bed_file_from_gff_file_prokka.py [your gff file from prokka]
```

For example:
```
make_bed_file_from_gff_file_prokka.py PROKKA_09222020.gff
```

This will create a bed file that ends in .bed. You can take a look at it if you wish-- it should have the contig name, the coordinates of your ORF, and the name of your ORF.



#### 18. Run script to calculate coverage
Now run a script that will use your bed file and will calculate the read depth for every single ORF in your ORF file. You might have to provide the path to your bed file, or copy it into your project directory.

```
samtools bedcov [your bed file] [your sorted bam file] > [ an output file that ends in _ORF_coverage.txt]
```

For example:
```
samtools bedcov PROKKA_09222020.bed ERR599166_mapped_sorted.bam > ERR599166_ORF_coverage.txt
```

The `samtools bedcov` command will give you a file that ends in `ORF_coverage.txt`.

#### 19. Matching the Prokka annotations with your ORF coverage

Now you have the coverage for all of your ORFs. To make it easier to read and work with the results, let's match the coverage of each ORF with the Prokka annotation of that ORF so you can more easily read and understand the results. To do that, run this script, substituting in the names and paths of your own ORF coverage file and your faa file from Prokka:
```
python /Accounts/Genomics_Bioinformatics_shared/python_scripts/merge_name.py [ORF_coverage.txt file] [ORFs.faa file]
```

So, for example:
```
python /Accounts/Genomics_Bioinformatics_shared/python_scripts/merge_name.py ERR598966_ORF_coverage.txt prokka_project/PROKKA_01252021.faa
```

Your output will be a txt file that matches the name of your ORF coverage file and ends in `_matched.txt.`

#### 20. Calculate coverage in Excel and examine results

Use FileZilla or `scp` to move the new `_matched.txt.` file to your local computer, then open with Excel. The output file should give the name of your contig, the start coordinate, the stop coordinate, the name of your open reading frame, the sum of the per-base coverage, and then the Prokka annotation for that protein. To get the *average* coverage for each ORF, divide the sum of the per-base coverage by the difference between the stop and start coordinates. In other words, type "average coverage" in the top of column H, then one row down, type this and then fill down to the bottom:
=F2/(D2-C2)

*Pro tip: Rather than drag for 16,000 rows in Excel, highlight the top cell, then scroll down to the bottom of the column, then hold 'Shift' while you click the bottom cell of the column where the data ends, then click Edit --> Fill--> Down.*

Now you have the average coverage of every ORF for this particular bam file, and you should be able to see what the annotation of each ORF is. You can explore your results and see what kinds of genes tend to have the highest and lowest coverage by sorting your spreadsheet in different ways.

This is a really common type of analysis for 'omics-based studies-- you can compare the coverage of specific genes of interest. For example, you might compare the coverage of genes related to photosythesis, respiration, and nitrogen fixation if you're interested in how abundant those metabolisms are in the community. We also use a very similar technique when we're doing expression analyses using something like RNASeq (more on that in coming weeks).

#### 21. Check for understanding

**Take a look at the spreadsheet and discuss with your lab group:**

Q5. As you scroll through the data file reporting the average coverage of all of your ORFs, which ORFs had the highest coverage? What do they encode? Speculate on why those genes may have had the highest coverage of all the genes in your dataset. NOTE! The genes with the highest coverage were probably really short and resulted from Illumina sequencing error-- i.e., ATATATATATATAT. IDBA-UD orders contigs by length, so I recommend skipping the contigs that are really short (high numbers in the contig name) and find the contig with the highest coverage that was unlikely to be a sequencing error.



**Another script that may be useful for your postlab or your final project, but isn't required for this lab:**

Let's say you ran a BLAST to identify ORFs of interest, and you want to only get the coverages of the genes that had a match in your BLAST results. To do that, run this script:
```
get_ORF_covg_from_BLAST_hits.py [BLAST file] [ORF coverage file]
```

#### 22. Copy your mapping files to the shared class folder

Please copy your bam files, bai files, bed files, and your ORF coverage files over to the class shared folder. Before doing that, you might want to change the names of some of your files so they are uniform and recognizable (see below example, and substitute your file names for the ones below).

For example:
```
mv PROKKA_09222020.bed ERR598995_assembly_ORFs.bed
```
Then:

```
cp ERR598995_mapped_sorted.bam /Accounts/Genomics_Bioinformatics_shared/mapping
cp ERR598995_mapped_sorted.bam.bai /Accounts/Genomics_Bioinformatics_shared/mapping
cp ERR598995_assembly_ORFs.bed /Accounts/Genomics_Bioinformatics_shared/mapping
cp ERR598995_ORF_coverage.txt  /Accounts/Genomics_Bioinformatics_shared/mapping
```


#### 23. This week's postlab writeup
For this week's post-lab writeup:

**Mini Research Question**

Write either a question or generate a hypothesis about the relative coverage of this set of genes with respect to your project datasets.


*Example #1:*
I hypothesize that there will be lower coverage of genes related to photosynthesis (i.e. the psb genes) in the mesopelagic zone relative to the surface. This is because at the surface there will be more organisms that photosynthesize compared to the mesopelagic zone, where less light is available. Therefore, a lower proportion of genes in the microbial community in the mesopelagic zone will be related to photosynthesis compared to the surface, and therefore, fewer reads will map to photosynthesis genes in the mesopelagic zone.

*Example #2:*
I hypothesize that there will be higher coverage of genes related to viruses in my Lyman Lakes sample relative to the Cannon River samples because there are more viruses in Lyman Lakes than in river water, simply because there are more organisms to infect in lake waters.


Once you've identified the set of genes related to a specific metabolism/function/type of organism, and you have written a question or generated a hypothesis, find the average coverage to each of those ORFs in your dataset. Remember that more than one ORF may have that function.

Many of you will probably want to compare your mapping of your own reads to your own dataset to a mapping made by one of your classmates to their own dataset. The bam files and ORF coverage files should be saved in /Accounts/Genomics_Bioinformatics_shared/mapping.

**Describe your results and create at least one graph to visualize those results. This should represent a mini 'Results' section in a lab report or paper. Interpret your results within the context of the ecosystem you are investigating. This should represent a mini 'Discussion' section in a lab report or paper.**

**Please submit on Moodle by lab time next week. As before, please put your student ID number rather than your name.**
