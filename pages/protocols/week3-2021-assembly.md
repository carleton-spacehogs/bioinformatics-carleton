# Week 3: Introduction to Unix/Linux and the Server; Assembly with IDBA-UD; ORF Calling and Annotation with Prokka

Rika Anderson,
Carleton College

## Connecting to baross

#### 1. About baross
We are going to do most of our computational work on a remote server (or computer) called baross, which is a remote server with 96 CPUs, 50 TB of storage, and 768 GB RAM. (In case anybody cares, baross is named after [this](https://depts.washington.edu/astrobio/wordpress/profile/john-baross/) absolute legend of a scientist, who happened to be my PhD advisor.) baross (the server) lives in the basement of the CMC. You can access it from lab computers on campus, and you can also access it from your own computer at home. First, you have to learn how to access baross.

**IF YOU ARE IN A COMPUTER LAB ON THE CARLETON CAMPUS:**
Boot as a Mac user on the lab computer.

**IF YOU ARE WORKING ON A PERSONAL MAC COMPUTER:**
You should be able to follow the directions below.

**IF YOU ARE WORKING ON A PERSONAL COMPUTER THAT IS OPERATING WINDOWS OR ANOTHER OS:**
You will need to find a way to connect to a remote server. I recommend installing a [Ubuntu terminal.](https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)
If that doesn't work, you can use [PuTTY.](https://www.howtogeek.com/311287/how-to-connect-to-an-ssh-server-from-windows-macos-or-linux/)

#### 2. Opening Terminal
If you're on a Mac, find and open the Terminal application (it should be in "Applications" in the folder called "Utilities"). If you're on a PC or other OS, open the window you'll use to ssh into a remote server (like Ubuntu or PuTTY or something similar).

*The terminal is the interface you will use to log in to the remote server for the rest of the course. It is also the interface we will be using to run almost all of the bioinformatics software we will be learning in this course. You can navigate through files on the server or on your computer, copy, move, and create files, and run programs using the Unix commands we learned earlier this week.*

#### 3. ssh
We will use something called ssh, or a secure socket shell, to remotely connect to another computer (in this case it will our class server, baross). Type the following (substituting *username* with your own Carleton username-- the same name as your email address):

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

Sometimes we will run processes that are so computationally demanding that they will run for several hours or overnight. Some processes can sometimes take several days, or even weeks. In order to run these processes, you want to be able to close out your session on the remote server without terminating your process. To do that, you use screen. Type: `screen -S test`

Your terminal will open a clean window. You are now in a screen session that you have called 'test'. Any processes you start while you are in a screen session will continue running even after you log off of the remote server.

Let's say you've typed the commands for a process, and now you're ready to exit your screen session to let it run while you do other things. To leave the screen session, type: `Control+A d`

This will "detach" you from the screen session and allow you to do other things. Your screen session is still running in the background.

To resume your screen session to check on things, type: `screen -r test`. (Now you've re-entered your screen session called 'test.')

To kill your screen session (this is important to tidy up and keep things neat when your process is finished!) type: `Control+A k`.
- The computer will say: `Really kill this window [y/n]`. You type: `y`.
- If this doesn't work, you can also type `exit` to terminate the screen session.


## Creating an assembly and evaluating it (toy dataset)

#### 9. mkdir

Let’s say we’ve gone out in to the field and taken our samples, we’ve extracted the DNA, and we’ve sent them to a sequencer. Then we get back the raw sequencing reads. One of the first things we have to do is assemble them. To do that, we’re going to use a software package called [IDBA-UD](https://github.com/loneknightpy/idba). Bioinformaticians love to debate about which assembler is best, but ultimately it usually depends on the nature of your own dataset. If you find yourself with data like this someday and you want to know which assembler to use, my advice is to try a bunch of them and then compare the results using something like Quast, which we’ll test below. For now, we’re going to use IDBA-UD, which I’ve found to be a good general-purpose assembler.

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

We will talk more about how assembly works in class. Briefly, De Bruijn graph assembly works by searching reads for overlapping words, or "kmers." IDBA-UD starts by searching for small kmers in the short reads to assemble those short reads into longer contigs. Then it uses those constructed contigs for a new round of assembly with longer kmers. Then it does the same with even longer kmers. It iterates through kmers of length 20, 40, 60, 80, and 100. Each "contig-" file gives the results of each iteration for a different kmer size. The "scaffolds.fa" file provides the final scaffolded assembly, which pulls contigs together using the paired-end data.

#### 15. Examine output

Let's examine the assembly output files. First, take a look at your final scaffold output:

`less scaffold.fa`

#### 16. Fasta file

You should see a fasta file with some very long sequences in it. When you're done looking at it, type: `q`.

## Create and evaluate assembly (project dataset)

We have a few extra steps to better understand the output of your toy dataset assembly (see below), but first we're going to get started on the assembly on your project dataset because it's going to take a while. Before we start working on your project datasets, I encourage you all to start a computational lab notebook. I recommend that you open a text document in BBEdit or a similar text editor and save it on the cloud somewhere. In that notebook, it will be crucial that you keep a record of all the commands you run, especially for your project datasets. That way, if you get a funny result, or if you want to repeat an analysis, you can always go back to your lab notebook to see exactly what command you ran. In my own computational lab notebook, I record the dates and I include notes to myself about the analysis, or observations of the graphs I make. I recommend using a text document (like in BBEdit) rather than Word or Google Docs because you can copy-paste Unix commands from a text editor without worrying about formatting. These notebooks are just for you-- I won't be checking them-- but now that we're about to get into the real analysis, you should start one for your project datasets and maintain them throughout the rest of this course.

#### 17. screen

Since the assembly will take about 20-30 minutes (possibly more), I recommend running this on screen. Type:

```bash
screen -S assembly
```

#### 18. Make and change directories

Make a new directory in your home folder called "project_directory," and then change directory into your newly created directory:

```bash
cd ~
mkdir project_directory
cd project_directory
```

#### 19. Copy project dataset

Next, copy your project dataset in to that directory. Your assigned project dataset is listed in an Excel file that is on Moodle.

For example, if you have dataset `ERR590988` from the Arabian Sea, you would do this:

```bash
cp /usr/local/data/Tara_datasets/Arabian_Sea/ERR598966_sample.fasta .
```

#### 20. Start assembly

Next, start your assembly. **Substitute the name of your own sample dataset for the dataset shown here.**

```bash
idba_ud -r ERR598966_sample.fasta -o ERR598966_assembly
```

#### 21. Detach from screen
After you've started the process, detach from screen:

```bash
Ctrl+A d
```

These assemblies will take some time: possibly 30 minutes or more. To check on the progress of your assemblies, you can either use `top` to see if your assembly is still running, or you can re-attach to your screen using the instructions above. In the meantime, you can return to evaluating your toy dataset assemblies (step 22 below).



#### 22. Evaluate assembly quality

All right, now that we've got that going, let's return to your toy datasets. We’re going to evaluate the quality of our assembly. One measure of assembly quality is **N50**, which is the contig length for which the collection of all contigs of that length or longer contains at least 50% of the sum of the lengths of all contigs. To do that, we will use an assembly evaluator called [Quast](http://bioinf.spbau.ru/quast). Run it on your toy dataset. (Note that the command is kind of long, so scroll right in the little green box below.)

```bash
cd ~/toy_dataset_directory/toy_assembly
quast.py contig-20.fa contig-40.fa contig-60.fa contig-80.fa contig-100.fa scaffold.fa
```

Here is what the commands mean:

1. Invoke the program `quast.py`
2. We tell it to run the program on all of the output files of your toy assembly. Every fasta file that we list after `quast.py` will be included in the analysis.

Quast called the output `quast_results` by default. You can change that using the `mv` command:
`mv quast_results toy_assembly_quast_evaluation`


#### 23. View output by copying to your local computer

Quast gives you some nice visual outputs, but it's easier to look at these if you copy them from the server to your local computer. For this, we will use `scp`, or 'secure copy.' It's a lot like `cp`, except that you're copying from a remote server to your own computer. We will use `scp -r` because we're copying a directory recursively, not just a file.

Open up another Terminal window on your own computer. Don't ssh into anything. Navigate to your home directory.

```bash
cd ~
```
Type this (substituting your own username below):

```bash
scp -r username@baross.its.carleton.edu:~/toy_dataset_directory/toy_assembly/toy_assembly_quast_evaluation .
```

That means: copy something from the remote server baross and put it right here on my local computer.

Find your newly copied folder on your own computer. Double-click on the file called “report.html” and discuss it with your group (below).

#### 24. Pause to check for understanding
Take a pause here and check in with the other folks at their table. Help them catch up if they aren't at this step yet. When you're all at this point, discuss the following questions. You'll submit these as part of your postlab writeup described at the end of this protocol document (write your group members' names on the writeup as well).

- a) How many reads did you start with in the raw data file for your toy dataset, prior to assembly? (Hint: refer to the Unix tutorial you did, and use grep to search for the '>' symbol that denotes each new sequence. The > symbol also means something in Unix, so be sure to use quotes around it!)
- b) Examine the Quast output. Describe and explain the pattern you observe in terms of N50 as the kmer size increases (from 'contig-20.fa' all the way to 'contig-100.fa').
- c) Examine the Quast output. How does the N50 of your scaffold file compare to your contig-100 file? Explain why.
- d) Examine the Quast output. Which assembly output file had the most contigs? The fewest? Explain why.


#### 25. Bookkeeping

For the sake of future bookkeeping, you’ll need to modify your scaffold assembly file so that another software program we will use in the future (called `anvi’o`) is happy. **You must do this for both your toy and project datasets when they are finished assembling!** (Spaces and other weird characters are so tricky sometimes.) Do this while you are inside the directory with your assembly files:

```bash
anvi-script-reformat-fasta scaffold.fa -o toy_dataset_assembled_reformatted.fa -l 0 --simplify-names
```


## Searching for and annotating open reading frames

**We'll be doing this for both the toy and project datasets, so you can move on to this even if your project dataset is still assembling.**

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

#### 26. Run Prokka
Now run Prokka on your toy assembly, which is located in the toy_assembly folder:
```
prokka ../toy_assembly/toy_dataset_assembled_reformatted.fa --outdir prokka_toy
```
This means you're invoking prokka on your reformatted toy dataset assembly, and you're putting it in a new directory called `prokka_toy.`

#### 27. View output in FASTA format
You should see a directory called `prokka_toy`. Use `cd` to go into that folder, then use the program `less` to look at `PROKKA_01252021.faa` (or something like that-- adjust according the date). You should see a fasta file with amino acid sequences. Each amino acid sequence is an open reading frame (ORF), or a putative gene that has been identified from your assembly.

The cool thing is that Prokka has also annotated your genes with their putative function. You can see that in each sequence title, which provides the name of the sequence assigned by Prokka (e.g. KGLPOPID_00002) and the putative function (e.g. Proine/betaine transporter). A lot of them will say `hypothetical protein`, which simply means that Prokka couldn't find a good match for that protein in public databases.

Note that the file that ends in `.ffn` contains the same thing, except the ORF sequences are in nucleotide format, not amino acid. And the file that ends in `.gbk` is one commonly used in the National Institute of Health (NIH) National Centers for Biotechnology Information (NCBI) database, which we'll be using a lot over the next few weeks.

#### 28. View output in tab-separated column (tsv) format
Let's look at one last output format, which is in tab-separated columns, and therefore best visualized in a spreadsheet application like Excel. Use `scp` to copy this file from the server to your own computer. Remember, type this into a Terminal window that is NOT logged on to the server.

```
scp username@baross.its.carleton.edu:~/toy_dataset_directory/ORF_finding/prokka_toy/PROKKA_09252020.tsv .
```

#### 29. Open tsv file
Find your file on your local computer, and open your `.tsv` file in Excel (or Google Sheets). The column headers are in columns as follows, left to right:
1. `locus_tag`: the name that Prokka assigned the ORF
2. `ftype` (feature type): whether the thing it found is an ORF (CDS) or a tRNA or rRNA gene or something else.
3. `length_bp`: length of the gene
4. `gene`: if there was a match to an identified gene, it will give the name of that gene.
5. `EC_number`: the Enzyme Commission assigns numbers to known genes according to their function. If there was a match to a known gene with an EC number, that is given here.
6. `COG`: the Clusters of Orthologous Groups (COG) database groups proteins into categories of known function. If your ORF matched a specific COG, that is listed here. There is a website describing all of the COGS [here](ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/listCOGs.html).
7. `product`: The putative annotation for your ORF.

This document can be very useful: now you can see ALL of the ORFs in your metagenome and their putative annotation! This file can be handy when you're working on your postlab and/or when you're working on your final project.


## Annotate your own project datasets
Now that you have the general lay of the land, you're now going to annotate your own project dataset assemblies using `Prokka`.

#### 30. Check on your project assemblies
Hopefully your project assemblies are done by now. Log back into your screen session (`screen -r assembly`). When IDBA-UD finishes, it will say:

```
invalid insert distance
Segmentation fault
```

This is because we gave IDBA-UD single-end reads instead of paired-end reads. Never fear! Despite appearances, your assembly actually did finish. However, the use of single-end reads instead of paired-end reads means that your assemblies will NOT have scaffolds, because we couldn’t take advantage of the paired-end data to make scaffolds. So, your final assembly is just the final contig file from the longest kmer, called `contig-100.fa.`

#### 31. Rename assembly

When your project assembly is done, change the name of your final assembly from `contig-100.fa` to a name that starts with your dataset number, followed by `_assembly.fasta`. For example, you might call it something like `ERR598966_assembly.fasta`. This way we will have standardized names and save ourselves future heartache. While in your project_directory, type the following (substitute the name of your own project dataset):
```
 mv ERR598966_assembly/contig-100.fa ERR598966_assembly.fasta
 ```
 (You learned how to do this already-- see your Unix cheat sheet.)


#### 32. Run anvi-script-reformat-fasta

Run `anvi-script-reformat-fasta` on your completed project assemblies (see step 25 above). For example:
```
anvi-script-reformat-fasta ERR598966_assembly.fasta -o ERR598966_assembly_reformatted.fa -l 0 --simplify-names
```


#### 33. Run Prokka on your project dataset

If you aren't in your project directory right now, move into it and start `prokka`. You're probably still in a screen session, so you can just go ahead and run Prokka. Run the following, again, substituting your own project dataset name.

```
prokka ERR598966_assembly_reformatted.fa --outdir prokka_project
Ctrl+A d
```


## Playing with "big data"
Congratulations! You now have annotations for your entire project dataset. This means you have annotations for *thousands* of genes in your dataset. Feeling overwhelmed? Feeling excited about the endless possibilities? Feeling hungry for some questionable yet tasty LDC food? If it's one of the first two-- welcome to big data! If it's the last one, maybe lab has gone on for too long.

We're almost done. But first we're going to see an example of how we might analyze this kind of data.

#### 34. COG categories
First, there's a text file on the server that assigns every single COG in the COG database to a specific gene category, like `Translation, ribosomal structure, and biogenesis`. You can see it by typing this:

```
less /usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt
```

#### 35. Getting COG categories of your genes
Let's say you want to know how many of the genes in your dataset are responsible for translation, how many are for energy metabolism, how many are viruses, etc. Fortunately, there is a magic script available for you to do that on the server. Change directories into wherever your Prokka results are and type this:

```
get_ORF_COG_categories.py [name of your Prokka tsv file]
```
Et voilà! You'll get a new file that ends in `_cog_categories.txt` with a list of all the different COG categories and the number of genes that fall into that category. Some of the COGs fall into multiple categories, denoted with, for example, `Signal_transduction_mechanisms/Transcription`. So now, you have a detailed list of all of the ORFs and their functions (the original Prokka .tsv file) AND a list of what general categories your ORFs fall into (the new 'cog_categories.txt' file). Imagine the possibilities!


#### 36. Share your data

Please copy the data you generated today and put them into the folder that we'll be sharing as a class, substituting the placeholders below with the name of your own file. This will be the shared folder for all of the data that you generate as a class.

Copy your assemblies to the shared class assemblies folder:
```
cp [name of your formatted, assembled file] /Accounts/Genomics_Bioinformatics_shared/assemblies
```

Copy your Prokka tsv and faa files to the shared class Prokka folder, and rename it with the name of your project dataset so it's recognizable. for example:
```
cp PROKKA_09222020.tsv /Accounts/Genomics_Bioinformatics_shared/PROKKA_results/ERR598995_ORFs.tsv
cp PROKKA_09222020.faa /Accounts/Genomics_Bioinformatics_shared/PROKKA_results/ERR598995_ORFs.faa
```

Copy your COG categories file into the shared class Prokka folder, and change the name so it's easily recognizable. For example:
```
cp PROKKA_09222020_cog_categories.tsv /Accounts/Genomics_Bioinformatics_shared/PROKKA_results/ERR590988_cog_categories.txt
```

Obviously, please substitute the names of your own files for the ones listed above.


#### 37. Exit

If you haven't killed your screen yet, you should do so. As you did above, while in your screen session, type:
`Ctrl A` and then `k`. The the computer will say: `Really kill this window [y/n]`. You type: `y`.

When you are all done with your work on the server, you can log off the server by typing:

```bash
exit
```


## Lab assignment this week

Write a mini research question based on your ORF annotations. If you like, you can compare your sample to someone else's. This doesn't have to be overly in-depth: ask a simple question that can be answered by the files you and/or your classmates have generated today.

For example, you might ask, "Is there a difference in the number of genes related to the mobilome/prophage in the surface ocean compared to the mesopelagic zone?"

Or for example, you might ask "What genes are present in the mesopelagic zones that are missing from the surface ocean?"

Or, "Are there more genes labeled as "hypothetical" genes in the Antarctic Ocean compared to the North Pacific?"

Or, "Within a single sample, are there more genes related to energy metabolism or replication/recombination/repair?"

Make a plot or table that illustrates the answer to your question. In a paragraph or so, explain your results and speculate as to why your results look they way they do, and what else you might want to investigate related to this idea in the future.

Submit this plus your "check for understanding" questions above (in a single document) on the class Moodle page by lab time next week. (Note: your "check for understanding" questions come from group discussion; your mini research question should be done independently, but you're free to talk to each other for ideas or help troubleshooting. The final product should represent your own work.)

 **I prefer to grade these blind, so please put your student ID number, rather than your name, on your assignment. (This applies to all future assignments as well.)**
