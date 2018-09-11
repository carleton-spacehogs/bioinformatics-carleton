# Week 6: Binning genomes with anvi'o

## Intro
This week in lab we’ll learn how to visualize your metagenomes and pull out individual genome bins. We aren’t going to use a toy dataset this week—we’re going straight into analysis with your metagenome datasets for your projects.

Genomes are disentangled from metagenomes by clustering reads together according to two properties: **coverage** and **tetranucleotide** frequency. Basically, if contigs have similar coverage patterns between datasets, they are clustered together; and if they have similar kmers appear over and over again, they will cluster together. When we cluster contigs together like this, we get a collection of contigs that are thought to represent a reconstruction of a genome from your metagenomic sample. We call these 'genome bins,' or 'metagenome-asembled genomes (MAGs).'

There is a lot of discussion in the field about which software packages are the best for making these genome bins. And of course, the one you choose will depend a lot on your dataset, what you’re trying to accomplish, and personal preference. I chose anvi’o because it is a nice visualization tool that builds in many handy features.

I am drawing a lot of information for this tutorial from the anvi'o website. If you'd like to learn more, see the link below.

http://merenlab.org/2016/06/22/anvio-tutorial-v2/

## Preparing your contigs database for anvi'o
#### 1. ssh tunnel
Boot your computer as a Mac and use the Terminal to ssh in to baross. (For now, you will have to run all anvi'o-related things on baross.) Anvio requires visualization through the server, so this week we have to create what is called an "ssh tunnel" to log into the server in a specific way. Substitute "username" below with your own Carleton login name.

NOTE: Each of you will be assigned a different port number (i.e. 8080, 8081, 8082, etc.). Substitute that number for the one shown below.

```
ssh -L 8080:localhost:8080 username@baross.its.carleton.edu
```

#### 2. Copying data
Make a new directory called “anvio” inside your project folder, then change into that directory.

Copy the following into your new directory:

1) assembled contigs with the cleaned up names (i.e. >c_000000000001)
2) all of your sorted .bam and .bai files that you already mapped to those contigs last week

```
mkdir project_directory/anvio
cd project_directory/anvio
cp [path to assembled, formatted contigs] .
cp [path to sorted .bam files] .
cp [path to sorted .bai files] .
```

#### 3. Get gene calls and annotations from Prokka
The first thing you have to do is make contigs database, which contains the sequences of your contigs, plus lots of information about those contigs. This includes information from Prokka-- so we have to take the information from Prokka and put it in our contigs database.

First, you hav to run a script on your Prokka files to convert them into a text file that we can import into anvi'o. Navigate to wherever your Prokka results are for your project assembly, and run this script:

```
gff_parser.py [your PROKKA gff file] --gene-calls gene_calls.txt --annotation gene_annot.txt
```
Copy the `gene_calls.txt` and `gene_annot.txt` files over to your anvi'o folder.

#### 4. Make the contigs database
Now, you make the database.

-`anvi-gen-contigs-database` is the anvi’o script that makes the contigs database.
-`–f` is the fasta file with your contigs that you have already assembled and fixed.
-`–o` provides the name of your new contigs database.
-`external_gene_calls` provides the name of the Prokka file you just made so you can import the Prokka calls into your contigs database
-`--ignore-internal-stop-codons` will ignore any internal stop codons in your gene calls. Sometimes these will get included in your Prokka results by accident, but for our purposes we can ignore them.

```
anvi-gen-contigs-database -f [your formatted, assembled contigs] -o contigs.db --external-gene-calls gene_calls_PROKKA.txt --ignore-internal-stop-codons
```

#### 5. Import the Prokka annotations
 Import the functional annotations like this:
 ```
 anvi-import-functions -c contigs.db -i gene_annot.txt
 ```

#### 6. Search for single copy universal genes
Now we will search our contigs for archaeal and bacterial single-copy core genes. This will be useful later on because when we try to disentangle genomes from this metagenome, these single-copy core genes can be good markers for how complete your disentangled genome is.

This process is slow, so we're going to run it on 5 CPUs rather than just 1. You can run it on screen in the background while you move forward with step 6. It should take a little under 10 minutes.

```
screen
anvi-run-hmms -c contigs.db -T 5
```

#### 7. Determine taxonomy using Centrifuge
Now we are going to figure out the taxonomy of our contigs using a program called centrifuge. Centrifuge is a program that compares your contigs to a sequence database in order to assign taxonomy to different sequences within your metagenome. We're going to use it first to classify your contigs.

If you would like to know more, go here: http://merenlab.org/2016/06/22/anvio-tutorial-v2/ and here: http://www.ccb.jhu.edu/software/centrifuge/

First, export your genes from anvi'o.
```
anvi-get-dna-sequences-for-gene-calls -c contigs.db -o gene-calls.fa
```

#### 8. Run Centrifuge
```
centrifuge -f -x /usr/local/CENTRIFUGE/p_compressed gene-calls.fa -S centrifuge_hits.tsv
```

#### 9. Import Centrifuge data
Now import those centrifuge results for your contigs back in to anvi'o. It has a parser written into the software that can automatically read and import centrifuge output.
```
anvi-import-taxonomy -c contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge
```

## Incorporating mapping data

#### 10. Copy mapping files
In order to make bins, anvi'o needs to compare mappings from different datasets. Today, we are going to compare all of the datasets that are at the same depth as yours: for example, if you are assigned to the surface layer, you should pull in all the other surface layer mappings to your own dataset. Check the Moodle for the data spreadsheet explaining which is which. You need both the **sorted .bam files** and the **.bai files**. Those are stored at ``/Accounts/Genomics_Bioinformatics_shared/Tara_mappings/``. Copy the ones you need over to the directory that you are in now.

```
cp /Accounts/Genomics_Bioinformatics_shared/Tara_mappings/[bam and/or bai files you want] .
```

#### 11. Import mapping files into anvi'o with anvi-profile
Now anvi’o needs to combine all of this information—your mapping, your contigs, your open reading frames, your taxonomy—together. To do this, use the anvi-profile script.

-`anvi-profile` is the name of the program that combines the info together
-The `–i` flag provides the name of the sorted bam file that you copied in the step above.
-The `-T` flag sets the number of CPUs. There are 11 of you, and 96 to spare. For now, let's set it to 5 so we don't blow up the server.
-The `-M` flag sets a minimum contig length. In a project for publication, you'd want to use at least 1000, because the clustering of contigs is dependent on calculating their tetranucleotide frequencies (searching for patterns of kmers). You need to have a long enough contig to calculate these frequences accurately. But for our purposes, let's use 500 so you can use as many contigs as possible.
```
anvi-profile -i [your sorted bam file] -c contigs.db -T 5 -M 500
```

#### 12. Merge them together with anvi-merge
Now merge all of these profiles together using a program called anvi-merge. You have to merge together files in directories that were created by the previous profiling step. The asterisk * is a wildcard that tells the computer, 'take all of the folders called 'PROFILE.db' from all of the directories and merge them together.'

We're also going to tell the computer not to bin these contigs automatically (called 'unsupervised' binning), we want to bin them by hand ('supervised' binning). So we use the --skip-concoct-binning flag.

This step will take a couple minutes.
```
anvi-merge */PROFILE.db -o SAMPLES-MERGED -c contigs.db --skip-concoct-binning
```
## Visualizing and making your bins

#### 13. anvi-interactive
Now the fun part with pretty pictures! Type this to open up the visualization of your contigs:
```
anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
````

#### 14. Visualize in browser
Now, open up a browser (Chrome works well) and type this into the browser window to stream the results directly from the server to your browser. NOTE that each of you will be assigned a different port number (i.e. 8080, 8081, 8082, etc); use the one you logged in with in the first step.

http://localhost:8080

Cool, eh?

Click 'Draw' to see your results! You should see something like this:
![anvio screenshot](../images/anvio.png)

*What you are looking at:

-the tree on the inside shows the clustering of your contigs. Contigs that were more similar according to k-mer frequency and coverage clustered together.

-the rings on the outside show your samples. Each ring is a different sample. (So if you had a deep chlorophyll max sample, the rings are all different deep chlorophyll max samples.) This is a visualization of the mapping that you did last week, but now we can see the mapping across the whole sample, and for all samples at once. There is one black line per contig. The taller the black line, the more mapping for that contig.

-the 'taxonomy' ring shows the centrifuge designation for the taxonomy of that particular contig.

-the 'GC content' ring shows the average percent of bases that were G or C as opposed to A or T for that contig.*

#### 15. Make bins
We will go over the process for making bins together in class.

Because your datasets are fairly small, your bins are also going to be very small. Your percent completeness will be very low. Try to identify ~3-5 bins according to patterns in the mapping of the datasets as well as the GC content.

When you are done making your bins, be sure to click on 'Store bin collection', give it a name ('my_bins' works), and then click on 'Generate a static summary page.' Click on the link it gives you. It will provide lots of information about your bins. In the boxes under the heading 'taxonomy,' you can click on the box to get a percentage rundown of how many contigs in your bin matched specific taxa according to centrifuge, if any matched.

**Once you have completed your binning process, take a screenshot of your anvi'o visualization and save it as 'Figure 1.' Write a figure caption explaining what your project dataset is, and which datasets you mapped to your sample.**

#### 16. Finding bin information
You will find your new bin FASTA files in the directory called '~/project_directory/anvio/SAMPLES-MERGED/SUMMARY_my_bins'.

`bins_summary.txt` provides just that, with information about the taxonomy, total length, number of contigs, N50, GC content, percent complete, and percent redundancy of each of your bins. This is reflected in the summary html page you generated earlier when you clicked 'Generate a static summary page.'

If you go to the directory `bin_by_bin`, you will find a series of directories, one for each bin you made. Inside each directory is a wealth of information about each bin. This includes (among other things):

-a FASTA file containing all of the contigs that comprise your bin (i.e. `Bin_1-contigs.fa`)

-mean coverage of each bin across all of your samples (i.e. `Bin_1-mean-coverage.txt`)

-files containing copies of single-copy, universal genes found in your contigs (i.e. `Bin_1-Rinke_et_al_hmm-sequences.txt` and `Bin_1-Campbell_et_all-hmm-sequences.txt`

-information about single nucleotide variability in your bins-- the number of SNVs per kilobasepair. (i.e. `Bin_1-variability.txt`)

## Analyzing your bins
For this week's post-lab assignment, you won't do a mini-research question because the quality of your bins may vary, through no fault of your own. (That's real data for you...) So instead, we're going to answer a few questions about the data that you just generated.


##### 1. **Bin completeness**

The easiest bins to generate are often the ones that had the longest N50 and the lowest population diversity. Take a look at this file:

`SAMPLES-MERGED/SUMMARY_my_bins/bins_summary.txt`.

1a. Which bin was the most complete?

1b. What was the N50 of that bin?

1c. What was the anvi'o-labeled taxonomy of this bin?

1d. Now find your bin folder in `SAMPLES-MERGED/SUMMARY_my_bins/bin_by_bin`. Inside that folder is a file called `SAMPLES-MERGED/SUMMARY_my_bins/bin_by_bin/Bin_X/Bin_X-Campbell_et_al-hmm-sequences.txt`. Use `less` to open that file. These are the single-copy universal genes (like ribosomal genes, for example) that anvi'o used to characterize the completeness of this bin. As a bonus, we can use these sequences to help figure out what this bin is. Copy the longest sequence in that file and BLAST it using blastn on the NCBI webpage. Report the e-value and description of the top hit. Does it match the taxonomy that anvi'o reports?




##### 2. **Coverage across all bins for one sample**

Use the `SAMPLES-MERGED/SUMMARY_my_bins/bins_across_samples/mean_coverage.txt` file to answer these questions.

2a. Which bin had the highest coverage across all samples, when looking only at the mapping of your own dataset (self-to-self)?

2b. What is the taxonomy of that bin, according to anvi'o and according to BLAST (see part 1 above)?

2c. Make a bar graph in which each bin is a different bar, and the bar height indicates the mean coverage. Call this 'Figure 2' and include a figure caption.

2d. What does it mean, biologically/ecologically speaking, for a bin to have high coverage?



##### 3. **Coverage of one bin across samples**

Choose one bin that you're interested in. Use the file `SAMPLES-MERGED/SUMMARY_my_bins/bin_by_bin/Bin_X/Bin_X-mean_coverage.txt` to answer these questions.

3a. Make a bar graph in which each sample is a different bar, and the bar height indicates the coverage of your bin in that sample. Call this "Figure 3" and include a figure caption.

3b. In which sample does your bin have the highest coverage? Second highest? What might this imply about the abundance of this microbe in different ecosystems?




##### 4. **Variability of all bins for one sample**

Use the `SAMPLES-MERGED/SUMMARY_my_bins/bins_across_samples/variability.txt` file to answer these questions.

4a. Which bins had the highest and lowest variability (single nucleotide variants per kilobase pair (SNVs/kbp) across all samples, when looking only at the mapping of your own dataset (self-to-self)?

4b. What is the taxonomy of those bins, according to anvi'o and BLAST (see above)?

4c. What do you think the variability might imply about the microbial populations represented by those bins? What might it imply about the selection pressure on those populations? (We'll talk more about this in class on Wednesday, so you might want to wait until after then...)

     *Important caveat! If your sample had low coverage (i.e. less than 10), that may be skewing the variability results because if you don't have any reads mapping, there are no variants to count!*

**Compile Figures 1-3 and your answers to questions 1-4 and submit on the Moodle by lab time next week.**
