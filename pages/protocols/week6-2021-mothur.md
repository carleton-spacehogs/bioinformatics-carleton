# Week 6: Classifying taxonomy of short reads with mothur

## Introduction
Today we're going to learn how to use sequence data to assess diversity across samples by focusing on just the 16S rRNA gene. You'll recall that 16S rRNA is like a "barcode" gene that we can use to compare and classify who is there. It's often used in microbiome studies. Today we'll learn how to:
  - classify the distribution of different types (taxa) of microbes in different datasets (i.e. what % are *Procholorococcus*, what % are SAR11, what % are Thaumarchaeota, and so on)
  - calculate diversity metrics to determine which sample sites had higher richness or evenness.


Before you start, think about which datasets you'd like to compare to your own. Perhaps you want to compare a bunch of surface samples from different regions, or perhaps you want to compare the three depths from your region. This could also potentially feed right into your final project. Choose 3-5 samples to compare (including yours).

## Using mothur to profile the taxonomy of your dataset

#### 1. Log in
Once you have decided on what samples to compare, log on to baross using ssh on your Terminal. Then make a new directory and change directory into it.
```
mkdir ~/project_directory/taxonomy
cd project_directory/taxonomy
```

#### 2. Copy over data
The Tara Ocean people have already identified all of the reads that matched 16S ribosomal RNA from their metagenomes, and they made separate FASTA files with just those reads. (Thanks, French scientists!) I copied those FASTA files into the directories listed in `/usr/local/data/Tara_datasets/`. Copy over all of the relevant 16S rRNA files that you will need into your new taxonomy folder.
```
cp /usr/local/data/Tara_datasets/[project sample site of interest]/[16S rRNA data file of interest] .
```

For example, you might do something like this:
```
cp /usr/local/data/Tara_datasets/Arabian_Sea/ERR598966_MERGED_FASTQ_16SrRNA_10000.fasta .
cp /usr/local/data/Tara_datasets/coastal_South_Africa/ERR598972_MERGED_FASTQ_16SrRNA.fasta .
cp /usr/local/data/Tara_datasets/North_Pacific/ERR598995_MERGED_FASTQ_16SrRNA_10000.fasta .
```

***Hot tip!***
If you want to copy over all of the 16S rRNA datasets from one location, you could use the asterisk (wildcard). For example:
`cp /usr/local/data/Tara_datasets/North_Pacific/*16S* .`
That would copy over all of the files containing `16S` in their title in the North Pacific folder.


#### 3. Open mothur
Now, we're going to use a program called `mothur` to analyze these sequences. mothur is a bit different from other programs that we've used in that we can enter a mothur interface and type commands that are specific to mothur. To enter the mothur program, simply type this:
```
mothur
```

#### 4. Create groups file
Right now, the 16S rRNA-matching reads from each sample site are in separate files. Pretty soon we are going to merge all of your FASTA files together in order to compare them. Before we do that, we need to tell mothur how to tell them apart once they are merged. We do that with the make.group command. Replace the file names below with your 16S rRNA fasta file names. Substitute the group names with a suitable name for your sample so that you can recognize it, like `NPacific_DCM` or `Arabian_surface.`

*Note: any time you see something in [brackets] in a command, it means you have to substitute that with your own data. Don't include the brackets in your command.*

```
make.group(fasta=[file1.fasta]-[file2.fasta]-[file3.fasta], groups=[group1]-[group2]-[group3])
```

For example:
```
make.group(fasta=ERR598944_MERGED_FASTQ_16SrRNA_10000.fasta-ERR599001_MERGED_FASTQ_16SrRNA_10000.fasta-ERR599078_MERGED_FASTQ_16SrRNA_10000.fasta, groups=meso-transect-surface)
```

***Hot tip #2!*** You can use the command system() if you want to use Unix commands while you are using mothur. If you can't remember the names of the 16S files you just copied, you can see them by typing this:
```
system(ls)
```


#### 5. Look at the groups file
This command should generate a file that is called either `groups` or `merge.groups`. Take a look at it with `less`.

```
system(less groups)
```
You will see that each sequence name is linked up with the group name that you provided. That way mothur can combine all of the sequences together into one file, but you can still keep track of which one belongs to which sample. This file will be essential for allowing mothur to compare your samples later on.


### 6. Merge FASTA files together
Now you can merge all of your FASTA files together, and the .groups file will record which sequences came from which file. The output here will be a file called `merged.fa`. Again, substitute "file-1.fa" and so on with the names of your 16S rRNA fasta files.

***Hot tip #3!*** use the up arrow on your keyboard to call up the last command you typed, and then edit that command instead of retyping all of the filenames again.

```
merge.files(input=[file1.fa]-[file2.fa]-[file3.fa], output=merged.fa)
```

For example:
```
merge.files(input=ERR598944_MERGED_FASTQ_16SrRNA_10000.fasta-ERR599001_MERGED_FASTQ_16SrRNA_10000.fasta-ERR599078_MERGED_FASTQ_16SrRNA_10000.fasta, output=merged.fa)
```

#### 7. Classify your sequences
Everything we've done so far is basically record-keeping. Now it's time to do SCIENCE!

We will classify our sequences by comparing them to a reference database. We will the SILVA database to compare these sequences. (It's a very good, well-curated 16S rRNA database.) 

Note!
You will see lots of warnings along the lines of: "[WARNING]: xxx could not be classified." We are going to have to leave these sequences out of the analysis! This means all of the unknowns will be grouped together even though they most likely represent many different species, so they will be missing from our diversity analyses later on. (In class, we'll talk about why we would ideally use something like operational taxonomic units, or OTUs, to do this analysis. Unfortunately, the Tara metagenomic data is too messy to be able to make nice OTUs, so we're going to classify every sequence.)

```
classify.seqs(fasta=merged.fa, group=groups, reference=/usr/local/data/silva_databases/silva.seed_v119.align, taxonomy=/usr/local/data/silva_databases/silva.seed_v119.tax)
```

#### 8. Open classified sequences
Use scp to transfer your files over to your local desktop. Open the file that is called `merged.seed_v119.wang.tax.summary` in Excel. (You may have to change the name so the file ends in '.txt' or Excel won't recognize it as a valid file to open.) You have seen this dataset before-- it's the one we worked on with Lin!

Here is the definition of the columns, from left to right:

- Taxonomic level is in the farthest left-hand column. The lower the number, the larger the phylogenetic classification, starting with domain, then phylum, class, order, family, genus, species. For example, Archaea, Bacteria, Eukarya, and 'unknown' are taxonomic level 1. Taxonomic level 2 classifies different phyla of Archaea, Bacteria, and Eukaryotes. Taxonomic level 3 classifies different classes of those phyla, and so on.
- The rankID provides a means of keeping track of where that particular organism falls. For example, the SAR_11 clade is rankID 0.2.17.2.9, which means it is a clade within the Alphaproteobacteria (0.2.17.2), which are a clade within the Proteobacteria (rankID 0.2.17), which is a clade within the Bacteria (rankID 0.2).
- The taxon column tells you the name of the taxon.
- The daughter level tells you how may levels down you are in the phylogeny.
- The 'total' tells you how many total sequences are within that taxonomic category.
- Each of the following columns gives you the taxonomic breakdown for that sample.

Part of your postlab assignment will be to explore this dataset and ask a scientific question about it.


#### 9. Share your data
Your classmates may wish to use your taxonomy data for their project datasets. Please rename your taxonomy files and share them on the class_shared directory. 

```
mv merged.seed_v119.wang.tax.summary [newname]
cp [newname] /Accounts/Genomics_Bioinformatics_shared/taxonomy
```
For example:
```
mv merged.seed_v119.wang.tax.summary rikas_taxonomy.summary
cp rikas_taxonomy.summary /Accounts/Genomics_Bioinformatics_shared/taxonomy
```

#### 10. Calculate the diversity of your sample
Now we're going to calculate the diversity at each of your sample sites. These will be pure calculations in Excel rather than on the command-line. mothur can do this using OTUs, but our data is too messy to create nice OTUs!

We'll calculate diversity using two measures: species richness (r) and the Shannon-Weiner Index (H'). The Shannon-Weiner Index (H') is meant to take into account both the taxon richness (i.e. how many different taxa there are) and evenness (i.e. does one taxon dominate, or are they evenly distributed?) 

We are going to calculate the diversity of each of your sample sites at taxonomic level 2. This means that each entry at level 2 (i.e. 'Euryarchaeota,' 'Thaumarchaeota,' 'Proteobacteria') will count as one taxon. Calculate the **richness** and the **Shannon-Weiner index** for each of your samples.

Richness = r = number of taxa in your sample

The equation for the Shannon-Weiner index is:

H’ = -Σ (Pi ln(Pi))

H’ = index of taxonomic diversity, the Shannon-Weiner Index

Pi = proportion (percent) of total sample belonging to the ith taxon

ln = natural log (log base e = not the same as log!)

This index takes into account not just the number of taxa (richness) in a sample, but also how evenly distributed the taxa are (evenness) within the sample.  The index increases by having more richness and/or by having greater evenness.

***Hint***: -Σ (Pi ln(Pi)) = -((Ptaxon1\*ln(Ptaxon1)) + (Ptaxon2\*ln(Ptaxon2)) + (Ptaxon3\*ln(Ptaxon3)) + …)

I suggest that you start by calculating the total number of sequences for each sample site for all of taxonomic level 2. Then, for each taxon within taxlevel2, calculate the total proportion of sequences belonging to that taxon (Pi). Then calculate H'.

Excel command for natural log: LN()
Note that in Excel, LN(0) = error, so skip the cells with a value of 0.
Pay attention to your use of parentheses!

Please do these calculations in a way that is clear so that I can track your calculations. You will be submitting these Excel spreadsheets as part of your lab assignment for this week. Please compile the total number of sequences, the species richness, and the Shannon-Weiner index for each of your sample depths in a table. **Save this as Table 1 and provide a caption.**

#### 11. Check for understanding

**Check in with your lab group. Help each other catch up, and then discuss the following check for understanding questions and write down your responses for your postlab:**

Q1. Many of your sequences were unclassifiable. How would this likely affect your richness calculations for each sample? Explain why.

Q2. What is the difference between richness and the Shannon-Weiner index? Describe a situation in which you might have a high richness but a relatively low Shannon-Weiner index.

#### 12. Share data to Google Sheet

Last, please enter the richness and Shannon-Weiner Index data for your sample site in the Google Sheet [here](https://docs.google.com/spreadsheets/d/1BOEHLLBYoymVR-7bITEpetpqa8yS_n2HNKojKAZLIQc/edit?usp=sharing) . We will use this data on Friday with Lin, and you may use it for your postlab this week.

#### 13. Mini research question
Design a mini research question using the data you've generated today. Remember that you have metadata (i.e. temp, chlorophyll, etc.) available as well in the Google Sheet. Generate a plot or set of plots addressing your research question. Each plot should have a figure caption. (On Friday, we will work on data visualization on these types of datasets with Lin.) Then write a couple paragraphs describing your results, as if this were a mini Results and Discussion section. What kinds of trends do you see? What were you expecting to see, and do your results support your initial expectations? Explain your results in light of what we have discussed in class (and perhaps based on what you have seen in your previous postlabs) about trends in taxonomy and/or diversity in various ocean basins and/or at different depths.

 **Summary of what to turn in next week:**
 For this week's post-lab assignment, please submit the following (should all be in the same document):
 - Table 1
 - Responses to 2 "Check for understanding" questions
 - Mini research question
