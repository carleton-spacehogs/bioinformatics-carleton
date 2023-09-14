# Week 3: Introduction to Unix/Linux and the Server; Assembly with IDBA-UD; ORF Calling and Annotation with Prokka

Rika Anderson,
Carleton College

## Connecting to baross

#### 1. About baross
We are going to do most of our computational work on a remote server (or computer) called baross, which is a remote server with 96 CPUs, 50 TB of storage, and 768 GB RAM. (In case anybody cares, baross is named after [this](https://depts.washington.edu/astrobio/wordpress/profile/john-baross/) absolute legend of a scientist, who happened to be my PhD advisor.) baross (the server) lives in the basement of the CMC. You can access it from lab computers on campus, and you can also access it from your own computer at home. First, we have to learn how to access baross.

**IF YOU ARE IN A COMPUTER LAB ON THE CARLETON CAMPUS:**
Boot as a Mac user on the lab computer.

**IF YOU ARE WORKING ON A PERSONAL MAC COMPUTER:**
You should be able to follow the directions below.

**IF YOU ARE WORKING ON A PERSONAL COMPUTER THAT IS OPERATING WINDOWS OR ANOTHER OS:**
You will need to find a way to connect to a remote server. Some Windows machines have a native Terminal now. If you don't have one, I recommend installing a [Ubuntu terminal.](https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)
If that doesn't work, you can use [PuTTY](https://www.howtogeek.com/311287/how-to-connect-to-an-ssh-server-from-windows-macos-or-linux/). 

#### 2. Opening Terminal
If you're on a Mac, find and open the Terminal application (it should be in "Applications" in the folder called "Utilities"). If you're on a PC or other OS, open the window you'll use to ssh into a remote server (like Ubuntu or PuTTY or something similar).

*The terminal is the interface you will use to log in to the remote server for the rest of the course. It is also the interface we will be using to run almost all of the bioinformatics software we will be learning in this course. You can navigate through files on the server or on your computer, copy, move, and create files, and run programs using the Unix commands we learned earlier this week.*

#### 3. ssh
We will use something called ssh, or a secure socket shell, to remotely connect to another computer (in this case it will be our class server, baross). Type the following (substituting *username* with your own Carleton username-- the same name as your email address):

```bash
ssh username@baross.its.carleton.edu
```

#### 4. ssh (cont.)

You should see something like this: `[your username]@baross.its.carleton.edu's password:`

Type in your Carleton password. NOTE!! You will not see the cursor moving when you type your password. Never fear, it is still registering what you type. Be sure to use the correct capitalization and be wary of typos.

#### 5. pwd

Whenever you ssh in to baross, you end up in your own home directory. Each of you has your own home directory on baross. To see where you are, print your current working directory.

```bash
# How to print your working directory (find out where you are in the system)
pwd
```


#### 6. Your path

You should see something like this: `/Accounts/[your username]`

This tells you what your current location is within baross, also known as your **path**. But before we actually run any processes, we need to learn a few important things about using a remote shared server.


## Important things to know about when using a remote shared server

#### 7. top

One of the most important things to learn when using a remote shared server is how to use proper **server etiquette.** You don't want to be a computer hog and take up all of the available processing power. Nobody earns friends that way. Similarly, if it looks like someone else is running a computationally intensive process, you might want to wait until theirs is done before starting your own, because running two intensive processes at the same time might slow you both down.

To see what other processes are being run by other users of the same server, type: `top`

You should see something like this:

![top screenshot](../images/top_screenshot.png)

Here, we can see that randerson is the only user on baross, I am using a process called top, and it's taking up 0.3% of one central processing unit (CPU), and zero memory. Not much is going on in this case. (All the other things listed under 'root' are processes that the computer is running in the background.) If someone’s process says 100.0 or 98.00 under %CPU, it means they’re using almost one entire CPU. **baross has 96 CPUs and 768 gigabytes of RAM total.** This is also my research server, so please be courteous and leave some CPUs for my research students. A good rule of thumb is to try to leave at least 15 CPUs available at any given time. If we try to overload the computer, it means that everyone else's processes will have to be shared among the available CPUs, which will slow everyone down. For example, if it looks like 80 CPUs are already being used at a particular time, you might want to wait until later to run your own process. Please also keep an eye on the memory (RAM) usage (where it says %MEM in top) and ensure we're not using more than 80% of total RAM across all jobs.

To quit out of top, type: `q`

#### 8. screen

Sometimes we will run processes that are so computationally demanding that they will run for several hours or overnight. In order to run these processes, you want to be able to close out your session on the remote server without terminating your process. To do that, you use screen. Type: `screen -S test`

Your terminal will open a clean window. You are now in a screen session that you have called 'test'. Any processes you start while you are in a screen session will continue running even after you log off of the remote server.

Let's say you've typed the commands for a process, and now you're ready to exit your screen session to let it run while you do other things. To leave the screen session, type: `Control+A d`

This will "detach" you from the screen session and allow you to do other things. Your screen session is still running in the background.

To resume your screen session to check on things, type: `screen -r test`. (Now you've re-entered your screen session called 'test.')

To kill your screen session (this is important to tidy up and keep things neat when your process is finished!) type: `Control+A k`.
- The computer will say: `Really kill this window [y/n]`. You type: `y`.
- If this doesn't work, you can also type `exit` to terminate the screen session.


## Creating an assembly and evaluating it (toy dataset)

#### 9. mkdir

Let’s say we’ve taken our samples, we’ve extracted the DNA, and we’ve sent them to a sequencer. Then we get back the raw sequencing reads. One of the first things we have to do is assemble them. To do that, we’re going to use a software package called [IDBA-UD](https://github.com/loneknightpy/idba). Bioinformaticians love to debate about which assembler is best, but ultimately it usually depends on the nature of your own dataset. If you find yourself with data like this someday and you want to know which assembler to use, my advice is to try a bunch of them and then compare the results using something like Quast, which we’ll test below. For now, we’re going to use IDBA-UD, which I’ve found to be a decent general-purpose assembler.

Make a new directory in your home folder:

```bash
mkdir toy_dataset_directory
```

#### 10. cd

Change directory into your toy dataset directory:

```bash
cd toy_dataset_directory
```

#### 11. Copy toy dataset

Copy the toy dataset into your assembly directory. Don’t forget the period at the end! This means you are copying into the directory you are currently in.

```bash
cp /usr/local/data/toy_datasets/toy_dataset_reads.fasta .
```

#### 12. Run idba-ud on toy dataset

Run idba-ud on the toy dataset. Here is what these commands mean:

1. Invoke the program `idba-ud`
2. The `-r` gives it the "path" (directions) to the reads for your toy dataset. `../` means it is in the directory outside of the one you're in.
3. The `-o` flag tells the program that you want the output directory to be called “toy_assembly”

```bash
idba_ud -r toy_dataset_reads.fasta -o toy_assembly
```

#### 13. cd to output directory

When that’s done running, go into your new output directory and take a look at what’s inside:

```bash
cd toy_assembly
ls
```

#### 14. The output files

You should see several files, including these:

```
contig.fa
contig-20.fa
contig-40.fa
contig-60.fa
contig-80.fa
contig-100.fa
scaffold.fa
```

IDBA-UD starts by searching for small kmers in the short reads to assemble those short reads into longer contigs. Then it uses those constructed contigs for a new round of assembly with longer kmers. Then it does the same with even longer kmers. It iterates through kmers of length 20, 40, 60, 80, and 100. Each "contig-" file gives the results of each iteration for a different kmer size. The "scaffolds.fa" file provides the final scaffolded assembly, which pulls contigs together using the paired-end data.

#### 15. Examine output

Let's examine the assembly output files. First, take a look at your final scaffold output:

`less scaffold.fa`

#### 16. Fasta file

You should see a fasta file with some very long sequences in it. When you're done looking at it, type: `q`.


#### 17. Evaluate assembly quality

Now that you have an assembly, we’re going to evaluate its quality. To do that, we will use an assembly evaluator called [Quast](http://bioinf.spbau.ru/quast). Run it on your toy dataset. (Note that the command is kind of long, so scroll right in the little green box below.)

```bash
cd ~/toy_dataset_directory/toy_assembly
quast.py contig-20.fa contig-40.fa contig-60.fa contig-80.fa contig-100.fa scaffold.fa
```

Here is what the commands mean:

1. Invoke the program `quast.py`
2. We tell it to run the program on all of the output files of your toy assembly. Every fasta file that we list after `quast.py` will be included in the analysis.

Quast called the output `quast_results` by default. You can change that using the `mv` command:
`mv quast_results toy_assembly_quast_evaluation`


#### 18. View output by copying to your local computer

Quast gives you some nice visual outputs, but it's easier to look at these if you copy them from the server to your local computer. For this, we have a couple of options. The easiest way is to use FileZilla, which should already be on the Mac computers in the computer lab. To do that, open the FileZilla application on your local computer, then enter the following:

Host: sftp://baross.its.carleton.edu
Username: Carleton username
Password: Carleton password
Port: 22
Click QuickConnect

You'll see your local computer contents show up on one side, and the stuff from the server show up on the other side. On the server side, nigate to wherever you kept the file you want to copy over. On the local computer side, navigate to wherever you want to put it. Double-click on the file you want to copy over. It should transfer over automatically. You can also drag and drop if you prefer.

Alternatively, if you don't have something like FileZilla, you can use a command called  `scp`, or 'secure copy.' It's a lot like `cp`, except that you're copying from a remote server to your own computer. We will use `scp -r` because we're copying a directory recursively, not just a file.

Open up another Terminal window on your own computer. Don't ssh into anything. Navigate to your home directory.

```bash
cd ~
```
Type this (substituting your own username below):

```bash
scp -r username@baross.its.carleton.edu:~/toy_dataset_directory/toy_assembly/toy_assembly_quast_evaluation .
```

That means: copy something from the remote server baross and put it right here on my local computer.


Regardless of the method you use, copy over the `toy_assembly_quast_evaluation` folder to your local computer. Double-click on the file called “report.html” and discuss it with your group (below).

#### 19. Pause to check for understanding
Take a pause here and check in with the other folks at their table. Help them catch up if they aren't at this step yet. When you're all at this point, discuss the following questions. Please put them into a shared Google Doc with your group and share that Google Doc with Rika. Please put everyone's name on the doc. I will not be grading these per se, but I'll be checking that you did them. This is a way to make sure everyone is on the same page and understanding what's going on.

- a) Examine the Quast output. Describe and explain the pattern you observe in terms of N50 as the kmer size increases (from 'contig-20.fa' all the way to 'contig-100.fa').
- b) Examine the Quast output. How does the N50 of your scaffold file compare to your contig-100 file? Explain why.
- c) Examine the Quast output. Which assembly output file had the most contigs? The fewest? Explain why.



## Searching for and annotating open reading frames

Now that we have assembled our contigs, we want to find genes on them. That requires identifying open reading frames (ORFs), and then comparing the ORF sequences with existing databases to see if we can figure out what kinds of genes they are. There are lots of programs to identify and annotate ORFs. We're going to use a program called [Prokka](https://github.com/tseemann/prokka), which wraps a bunch of different software into a pipeline that is nice and streamlined. Prokka basically does three things:

1. Identify ORFs
2. Compare those ORFs to databases of known genes to find the closest match and assign their putative functions
3. Several other things (like finding rRNA or tRNA genes, CRISPRs, etc.) that we aren't going to worry about right now.

Prokka works well and it's free, so we're using it today. It works best for archaeal and bacterial genomes. If you want to identify and annotate ORFs for eukaryotic genomes, there's lots of similar software out there to do that.

Go to your toy dataset directory and make a new directory:

```
cd ~/toy_dataset_directory
mkdir ORF_finding
cd ORF_finding
```

#### 20. Run Prokka
Now run Prokka on your toy assembly, which is located in the toy_assembly folder:
```
prokka ../toy_assembly/scaffold.fa --outdir prokka_toy --prefix toy
```
This means you're invoking prokka on your toy dataset assembly, and you're putting it in a new directory called `prokka_toy.` The prefix sets it up so that all of the files will start with `toy'.

#### 21. View output in FASTA format
You should see a directory called `prokka_toy`. Use `cd` to go into that folder, then use the program `less` to look at `toy.faa` (or something like that-- adjust according the date). You should see a fasta file with amino acid sequences. Each amino acid sequence is an open reading frame (ORF), or a putative gene that has been identified from your assembly.

The cool thing is that Prokka has also annotated your genes with their putative function. You can see that in each sequence title, which provides the name of the sequence assigned by Prokka (e.g. KGLPOPID_00002) and the putative function (e.g. Proine/betaine transporter). A lot of them will say `hypothetical protein`, which simply means that Prokka couldn't find a good match for that protein in public databases.

Note that the file that ends in `.ffn` contains the same thing, except the ORF sequences are in nucleotide format, not amino acid. And the file that ends in `.gbk` is one commonly used in the National Institute of Health (NIH) National Centers for Biotechnology Information (NCBI) database, which we'll be using a lot over the next few weeks.

#### 22. View output in tab-separated column (tsv) format
Let's look at one last output format, which is in tab-separated columns, and therefore best visualized in a spreadsheet application like Excel. Use `scp` to copy this file from the server to your own computer. Remember, type this into a Terminal window that is NOT logged on to the server.

```
scp username@baross.its.carleton.edu:~/toy_dataset_directory/ORF_finding/prokka_toy/toy.tsv .
```

#### 23. Open tsv file
Find your file on your local computer, and open your `.tsv` file in Excel (or Google Sheets). The column headers are in columns as follows, left to right:
1. `locus_tag`: the name that Prokka assigned the ORF
2. `ftype` (feature type): whether the thing it found is an ORF (CDS) or a tRNA or rRNA gene or something else.
3. `length_bp`: length of the gene
4. `gene`: if there was a match to an identified gene, it will give the name of that gene.
5. `EC_number`: the Enzyme Commission assigns numbers to known genes according to their function. If there was a match to a known gene with an EC number, that is given here.
6. `COG`: the Clusters of Orthologous Groups (COG) database groups proteins into categories of known function. If your ORF matched a specific COG, that is listed here. There is a website describing all of the COGS [here](ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/listCOGs.html).
7. `product`: The putative annotation for your ORF.

This document can be very useful: now you can see ALL of the ORFs in your metagenome and their putative annotation! This file can be handy when you're working on your postlab and/or when you're working on your final project.



## Create and evaluate assembly (project dataset)

Cool! Now that you've learned how to assemble and annotate your toy datasets, now we're going to start working on your project datasets. Before we do that, I encourage you all to start a **computational lab notebook**. I like to use a text editor like BBEdit or something similar and save it on Google Drive or something. In that notebook, it will be crucial that you keep a record of all the commands you run for your project datasets. That way, if you get a funny result, or if you want to repeat an analysis, you can always go back to your lab notebook to see exactly what command you ran. In my own computational lab notebook, I record the dates, and I include notes to myself about the analysis, or observations of the graphs I make. I recommend using a text document (like in BBEdit) rather than Word or Google Docs because you can copy-paste Unix commands from a text editor without worrying about formatting. These notebooks are just for you-- I won't be checking them-- but now that we're about to get into the real analysis, you should start one for your project datasets and maintain them throughout the rest of this course.


#### 24. Make and change directories

Make a new directory in your home folder called "project_directory," and then change directory into your newly created directory:

```bash
cd ~
mkdir project_directory
cd project_directory
```

#### 25. Copy project dataset reads

Next, copy your project dataset reads into that directory. Your assigned project dataset is listed in a Google Sheets file that is on Moodle. Some of you will have Tara Oceans samples, and some of you will have samples from the Arb-- note that they're located in different directories. This should be indicated on the spreadsheet.

For example, if you have dataset `ERR590988` from the Arabian Sea, you would do this:

```bash
cp /workspace/data/Genomics_Bioinformatics_shared/Tara_Oceans/reads/Arabian_Sea/ERR598966_sample.fasta .
```

or if you have an Arb dataset from sample 1A, you'd do this:

```bash
cp /workspace/data/Genomics_Bioinformatics_shared/Arb_sequences/reads/Arb_seqs_BIOL338_2023_GB_Fall2022_1A_merged.fa_randomsubset.fasta .
```

#### 26. Start assembly (except not really)

OK, so normally, I'd tell you to start your assembly using IDBA-UD. You would run a command like this:

```bash
idba_ud -r ERR598966_sample.fasta -o ERR598966_assembly
```

BUT DON'T DO IT!! Instead, I've pre-assembled your project datasets for you, because otherwise you'd be sitting here for hours. Copy them into your directory like this:

```bash
cp /workspace/data/Genomics_Bioinformatics_shared/Tara_Oceans/Tara_assemblies/your-assembly-name .
```

OR

```bash
cp /workspace/data/Genomics_Bioinformatics_shared/Arb_sequences/assemblies/your-assembly-name .
```

Now you have some pre-made assemblies ready for you to use.


## Annotate your own project datasets
Now you're going to annotate your own project dataset assemblies using `Prokka`. 

#### 27. Run Prokka on your project dataset

Normally, I'd tell you to open a screen session and run Prokka on your project dataset, just like you did for your toy dataset. You would run something like this:

```
screen -S prokka
prokka ERR598966_assembly_reformatted.fa --outdir prokka_project --prefix North_Pacific_prokka
Ctrl+A d
```

BUT DON'T DO IT!! As before, I've pre-run Prokka for all of you so that you don't have to sit here for hours waiting for it to finish. As above, please copy the Prokka folder into your own directory. So for example, you can do:


```bash
cp -r /workspace/data/Genomics_Bioinformatics_shared/Tara_Oceans/Tara_assemblies/prokka/your-Prokka-name .
```

OR

```bash
cp -r /workspace/data/Genomics_Bioinformatics_shared/Arb_sequences/prokka/your-Prokka-name .
```


## Playing with "big data"
Congratulations! You now have annotations for your entire project dataset. This means you have annotations for *thousands* of genes in your dataset. Feeling overwhelmed? Feeling excited about the endless possibilities? Feeling hungry for some questionable yet tasty LDC food? If it's one of the first two-- welcome to big data! If it's the last one, maybe lab has gone on for too long.

We're almost done. But first we're going to see an example of how we might analyze this kind of data.

#### 28. COG categories
First, there's a text file on the server that assigns every single COG in the COG database to a specific gene category, like `Translation, ribosomal structure, and biogenesis`. You can see it by typing this:

```
less /usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt
```

#### 29. Getting COG categories of your genes
Let's say you want to know how many of the genes in your dataset are responsible for translation, how many are for energy metabolism, how many are viruses, etc. Fortunately, there is a magic script available for you to do that on the server. Change directories into wherever your Prokka results are and type this:

```
get_ORF_COG_categories.py name-of-your-Prokka-tsv-file
```
Et voilà! You'll get a new file that ends in `_cog_categories.txt` with a list of all the different COG categories and the number of genes that fall into that category. Some of the COGs fall into multiple categories, denoted with, for example, `Signal_transduction_mechanisms/Transcription`. So now, you have a detailed list of all of the ORFs and their functions (the original Prokka .tsv file) AND a list of what general categories your ORFs fall into (the new 'cog_categories.txt' file). Imagine the possibilities!


#### 30. Share your data

Please copy the data you generated today and put them into the folder that we'll be sharing as a class, substituting the placeholders below with the name of your own file. This will be the shared folder for all of the data that you generate as a class.

You don't need to copy your assemblies and Prokka folders, since those were copied to begin with.

Copy your COG categories file into the shared class Prokka folder, and change the name so it's easily recognizable. For example:
```
cp PROKKA_09222020_cog_categories.tsv /Accounts/Genomics_Bioinformatics_shared/PROKKA_results/ERR590988_cog_categories.txt
```

Obviously, please substitute the names of your own files for the ones listed above.


#### 31. Exit

If you haven't killed your screen yet, you should do so. As you did above, while in your screen session, type:
`Ctrl A` and then `k`. The the computer will say: `Really kill this window [y/n]`. You type: `y`.

When you are all done with your work on the server, you can log off the server by typing:

```bash
exit
```


## Lab assignment this week

Write a mini research question based on your ORF annotations. If you like, you can compare your sample to someone else's. This doesn't have to be overly in-depth: ask a simple question that can be answered by the files you and/or your classmates have generated today.

For example, you might ask, "Is there a difference in the percent of genes related to the mobilome/prophage in the surface ocean compared to the mesopelagic zone?"

Or for example, you might ask "What genes are present in Lyman Lakes that are missing from the Cannon River?"

Or, "Is there a higher percentage of genes labeled as "hypothetical" in the Antarctic Ocean compared to Lyman Lakes?"

Or, "Within a single sample, are there more genes related to energy metabolism or replication/recombination/repair?"

Make a plot or table that illustrates the answer to your question. In a paragraph or so, explain your results and speculate as to why your results look they way they do, and what else you might want to investigate related to this idea in the future.

Submit this on the class Moodle page by lab time next week. 

 **I prefer to grade these blind, so please put your student ID number, rather than your name, on your assignment. (This applies to all future assignments as well.)**