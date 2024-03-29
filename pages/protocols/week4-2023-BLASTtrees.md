# Week 4: Local alignments and sequence search with BLAST, global alignments with MUSCLE, and making trees with RAxML

Rika Anderson,
Carleton College


## Logging in to the remote server

#### 1. Login
Boot as a Mac on the lab computer.

#### 2. Terminal
Find and open the Terminal application (`Applications/Utilities`). (For future reference if you're doing this from home, if you're on a PC , open your Ubuntu terminal or your PuTTY terminal.)

#### 3. Connect to baross

Log on to the server using your Carleton username and password.

```bash
ssh [username]@baross.its.carleton.edu
```


## Using BLAST

Last week, we annotated thousands of ORFs in one go. As awesome as that was, you may have noticed that there were a lot of proteins labeled as "hypothetical." Let's see if we can learn more about some of those genes by using BLAST, which is a tool that every bioinformatician should have in their toolkit.

#### 4. Copy a mystery sequence

Will use BLAST to compare your sequences against the National Center for Biotechnology Information (NCBI) non-redundant database of all proteins. This is a repository where biologists are required to submit their sequence data when they publish a paper. It is very comprehensive and well-curated.

Here's a mystery gene. Let's BLAST it. First, copy this sequence below.

```
>mystery_gene
MVPQTETKAGAGFKAGVKDYRLTYYTPDYVVRDTDILAAFRMTPQLGVPPEECGAAVAAESSTGTWTTVW
TDGLTSLDRYKGRCYDIEPVPGEDNQYIAYVAYPIDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRI
PPAYVKTFVGPPHGIQVERDKLNKYGRGLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQP
FMRWRDRFLFVAEAIYKAQAETGEVKGHYLNATAGTCEEMMKRAVCAKELGVPIIMHDYLTGGFTANTSL
AIYCRDNGLLLHIHRAMHAVIDRQRNHGIHFRVLAKALRMSGGDHLHSGTVVGKLEGEREVTLGFVDLMR
DDYVEKDRSRGIYFTQDWCSMPGVMPVASGGIHVWHMPALVEIFGDDACLQFGGGTLGHPWGNAPGAAAN
RVALEACTQARNEGRDLAREGGDVIRSACKWSPELAAACEVWKEIKFEFDTIDKL
```

#### 5. Navigate to NCBI site

Navigate your web browser to the [BLAST suite home page](https://blast.ncbi.nlm.nih.gov/Blast.cgi). Select Protein BLAST (blastp).

#### 6. Paste sequence and BLAST

Paste your sequence in to the box at the top of the page, and then scroll to the bottom and click “BLAST.” Give it a few minutes to think about it.

#### 7. Run tblastn

While that's running, open a new tab in your browser, navigate to the BLAST suite home page, and try blasting your protein using `tblastn` instead of `blastp`.

*What’s the difference?*

`blastp` is a protein-protein blast. When you run it online like this, you are comparing your protein sequence against the National Centers for Biotechnology Information (NCBI) non-redundant protein database, which is a giant database of protein sequences that is “non-redundant”—that is, each protein should be represented only once.

In contrast, `tblastn` is a translated nucleotide blast. You are blasting your protein sequence against a translated nucleotide database. When you run it online like this, you are comparing your protein sequence against the NCBI non-redundant nucleotide database, which is a giant database of nucleotide sequences, which can include whole genomes.

Isn't this a BLAST?!?
(Bioinformatics jokes! Not funny.)

#### 8. Pause to check for understanding with your lab group

Take a pause here and check in with the others at your table. If others haven't quite caught up to you yet, help them catch up. As a table, discuss and answer the following questions. As before, start a group Google Doc to record your answers and share it with Rika.

Q1. First look at your blastp results. Take a look at the list of hits below the top hit. Why do you think there are so many different hits with a similar function?

Q2. Describe the difference in your `tblastn` and `blastp` results. Why do they look different? How does this reflect the databases you are BLASTing against? Discuss and explain in what scenarios you might choose to use one over the other.


### Doing a BLAST against your own dataset
We just learned how to compare a sequence of interest against all proteins in the NCBI database. But let's say you just isolated a gene sequence that you're interested in, and want to know if any of the ORFs or contigs in **your** dataset has a match to just that sequence. In that case, you would want to blast against your own sequence set, not against the giant NCBI database. Here's how you do that.

First, make sure you are in the correct directory.

```
cd ~/toy_dataset_directory/ORF_finding/prokka_toy
```

#### 9. Make blast database
Turn your ORF dataset into a BLAST database. Here is what the flags mean:

- `makeblastdb` invokes the command to make a BLAST database from your data
- `-in` defines the file that you wish to turn into a BLAST database
- `-dbtype`: choices are "prot" for protein and "nucl" for nucleotide.

```
makeblastdb -in toy.faa -dbtype prot
```

#### 10. Find gene
Your ORFs are ready to be BLASTed. Now we need a gene of interest to search for. Let's say you are interested in whether there are any DNA polymerase genes here. (Fun fact: the toy dataset includes sequences from a thermophilic archaea from Yellowstone, which is where the DNA polymerase for PCR comes from!) The KEGG database is a handy place to find sequences that represent your gene of interest. Navigate your browser to [this KEGG website](https://www.genome.jp/kegg/).

You're seeing the main website for the Kyoto Enyclopedia of Genes and Genomes (KEGG), and it will be your friend for the rest of this class. This is a massive database with thousands of genes, organized according to category as well as by metabolic pathway.

Click on "KEGG Pathway," then click "Genetic Information Processing", then under "Replication and Repair" click "DNA replication." (You should land on [this page](https://www.genome.jp/pathway/map03030).)

You'll see a beautiful diagram of a DNA replication fork and all of the enzymes involved in this process. Next to the diagram, you'll see rectangles with gene names in them. Click on "PolB" under "DNA Polymerase B," which is next to the archaeal version of the replication fork. (Remember, the toy dataset happens to be archaeal.)

You'll see a table with lots of information about these genes. You'll also see a list of example genes under the table row called "Genes." Click on the first one (XNE: XNC1_4058), which happens to be a version of this gene that is found in a microbe called Xenorhabdus nematophila.

Scroll down, and you'll see some amino acid (AA) and nucleotide (NT) sequences.

#### 11. Make a FASTA file
Copy the amino acid sequence (don't copy the name, just the sequence itself). Then in the Terminal type "nano DNApol.fasta" and hit Enter. Type a ">" symbol, the name of the sequence and then paste the sequence in. It should look something like this:

```
>DNA polymerase
MPTVRGAEQGSKKRYAGLSGDKVIFKGLETVRTDWTPLAQTFQKELYTLIFHQQPYQEYI
REYIAKTMAGEFDDRLVYRKQLRRKLSDYQRNVPPHVRAVRLADEYNQQHNRPLQYQNGG
WINYLMTLSGPQPLENQTAAPDYQHYINKQLMPIADAILPFIQDNFMTLQTGQMNMTFE
```

After you've done this, type Ctrl+O, then hit Enter or Return, then Ctrl+X. This will close the text editor, and now you have a file called `DNApol.fasta`!


#### 12. BLAST it
Now, BLAST it! There are many parameters you can set. The following command illustrates some common parameters.

```bash
blastp -query DNApol.fasta -db toy.faa -evalue 1e-05 -outfmt 6 -out DNApol_vs_toy_ORFs.blastp
```

- `blastp` invokes the program within the blast suite that you want. (other choices are `blastn`, `blastx`, `tblastn`, `tblastx`.)
- `-query` defines your blast query-- in this case, the Pfam seed sequences for the CRISPR RAMP proteins.
- `-db` defines your database-- in this case, the toy assembly ORFs.
- `-evalue` defines your maximum e-value, in this case 1x10-5
- `-outfmt` defines the output format of your blast results. option 6 is common; you can check out [this link](https://www.ncbi.nlm.nih.gov/books/NBK279675/) for other options.
- `-out` defines the name of your output file. I like to title mine with the general format `query_vs_db.blastp` or something similar.

As we discussed in class, the e-value is a way to establish a cutoff of what makes a "good" BLAST hit. Smaller e-value = better hit.


#### 13. Examine results
Now let's check out your blast results. Take a look at your output file using `less`. (For easier visualization, you can also either copy your results (if it’s a small file) and paste in Excel, or transfer the file using `scp` and open it in Excel.)


#### 14. Column descriptions
Each blast hit is listed in a separate line. The columns are tabulated as follows, from left to right:

1. query sequence name
2. database sequence name
3. percent identity
4. alignment length
5. number of mismatches
6. number of gaps
7. query start coordinates
8. query end coordinates
9. subject start coordinates
10. subject end coordinates
11. e-value
12. bitscore



#### 15. BLAST comparison
Let’s try this another way. Run your blast again, but this time use a bigger e-value cutoff.

```bash
blastp -query DNApol.fasta -db toy.faa -evalue 0.2 -outfmt 6 -out DNApol_vs_toy_ORFs_evalue0.2.blastp
```


#### 16. Check for understanding

As before, take a pause and check in with your lab group. Help them catch up. As a group, take a look at your BLAST results and answer the following questions on the Google Doc.

Q4. How do these BLAST results differ from your previous BLAST? Explain why.

## Making An Alignment

OK, so we've learned how to do sequence search using BLAST, which is a *local* alignment search tool. Now we're going to learn how to make global alignments using MUSCLE, and then use those *global* alignments to make phylogenetic trees in order to learn about how genes are related.

We’ll start by creating a multiple sequence alignment with MUSCLE, and then we will make bootstrapped maximum likelihood phylogenetic trees with RAxML. We’ll use the Newick files generated by RAxML to visualize trees in an online tree visualization tool called the Interactive Tree of Life (iToL). We’ll do this using toy datasets, and then try out some trees on your project datasets.

#### 17. Aligning sequences
First we have to make a multiple sequence alignment with the sequences we wish to make into a tree. This could include any gene of interest. Today we're going to align a toy dataset made of DNA polymerase genes. I downloaded this from a database that contained many DNA polymerase genes, and I happened to include the sequence that had a match from the BLAST you just did.

Make a new directory within your toy dataset directory for making alignments and trees, then copy the toy dataset from the data directory to your toy dataset directory.
```
mkdir toy_dataset_directory/alignments_and_trees
cd toy_dataset_directory/alignments_and_trees
cp /workspace/data/Genomics_Bioinformatics_shared/DNA_polymerase.faa .
```

Take a look at it with `less`. It's just a FASTA file with a bunch of sequences.

#### 18. Align with MUSCLE
Now, we make a multiple sequence alignment using `muscle`. `muscle` uses dynamic programming to make global alignments of multiple sequences, like we did in class. It’s remarkably easy and fast.
```
muscle -in DNA_polymerase.faa -out DNA_polymerase.afa
```

What this means:

`muscle` is the name of the alignment program

`-in` defines the name of your input file, which can be either DNA or protein

`-out` defines the name of your output file. I like to give them an easy-to-recognize name with the extension `.afa`, which stands for “aligned fasta file.”

#### 19. Copy to a local computer using FileZilla OR scp
Let’s take a look at your alignment. I like to use the program Seaview to do this, and Seaview should be on the lab computer. You can use FileZilla to do this, or you can copy your file to your local computer using `scp.` (As before, substitute `[username]` with your own username.)

 ```
 scp [username]@baross.its.carleton.edu:/Accounts/your-username/toy_dataset_directory/alignments_and_trees/DNA_polymerases_COG_aligned.afa ~/Desktop
```
This lets you securely copy the aligned protein file from baross to your local computer. 

If you wanted to securely copy any file from your local computer to baross, use the command below. It's easy to find the path of where you want to put things-- simply navigate to where you want to put it on the server, and then either copy everything before the $ or just type `pwd`.

```
scp ~/Desktop/some_file.txt your-username@baross.its.carleton.edu:/path_of/your_destination_directory
```


#### 20. Visualize with Seaview
Open the application called “Seaview” and drag your file over to the alignment window. You should see something like this.

![seaview screenshot](../images/seaview.png)

Seaview shows the names of the sequences to the left. The letters to the right are the amino acids in your sequence, color-coded to make the alignment easier to see. You can easily see that some regions of the sequence are more highly conserved than others, and that some species appear to have an insertion or deletion in specific regions of the sequence. Note that this is an amino acid alignment, not a nucleotide alignment. (You could easily use MUSCLE to align nucleotides as well.)


## Making a Tree
Now we’re going to turn this alignment into a phylogenetic tree. We’re going to use a software package called RAxML, which is a commonly used tree-building software package that uses the maximum likelihood method to build trees.

#### 21. To make your tree, type this:

```
raxmlHPC-PTHREADS-AVX -f a -# 20 -m PROTGAMMAAUTO -p 12345 -x 12345 -s DNA_polymerase_aligned.afa -n DNA_polymerase.tree -T 4
```

What this means:

`-raxmlHPC-PTHREADS-AVX` is the name of the software package. This one is configured for the processors that are specific to this server.

`-f a` allows for rapid bootstrapping. This means RAxML will do a maximum likelihood search based on your protein sequences, and then bootstrap it as many times as you wish.

`-# 20` tells the program to bootstrap 20 times.

`-m PROTGAMMAAUTO` tells the program that these are protein sequences, and tells the program how to model protein evolution. To get into this is beyond the scope of this class, but fortunately RAxML is able to automatically choose the best one for us based on the number of sequences and the type of data we have.

`-p` and `–x` provide seed numbers so that the program can generate random numbers for the bootstrapping process.

`-s` gives the name of your aligned FASTA file

`-n` gives the name of your output Newick file, which will be made into a tree.

`-T` determines the number of threads. This is sort of like determining how many processors you'll use up for this process. Today, we'll use 4. Please don't use more than this without asking first.


NOTE: "bootstrapping" is a statistical test used to assess the reliability of your tree topology (its shape). We'll talk about this more in class next week.

#### 22. Tree output

You’ve made a tree! Congratulations! Let’s look at the raw RAxML output. You should have some files called:

`RAxML_bestTree.DNA_polymerase.tree`
`RAxML_bipartitionsBranchLabels.DNA_polymerase.tree`
`RAxML_bipartitions.DNA_polymerase.tree` <-- this is the one you want, because it gives you bootstrap values at the nodes on the tree.
`RAxML_bootstrap.DNA_polymerase.tree`
`RAxML_info.DNA_polymerase.tree`

Take a look at the bipartitions file using `less`.

This is a Newick file, and it’s a common format for phylogenetic trees.

#### 23. Visualizing the tree
Copy your Newick file (RAxML_bipartitions.DNA_polymerases_COG_aligned.tree) to your local computer using FileZilla or `scp`. Open up a web browser on your local computer and navigate to the [IToL website](http://itol.embl.de/) and use the class account:

username: biol338carleton
password: schiller

Click on "My trees" at the top and then click the + sign on the right to make a new workspace. Then click on “Upload tree files” and find your tree file and upload it.

Click on your tree. You should see it open in your window. You can play around with the settings on the right—you can make it circular or normal, you can display the bootstrap values as symbols or as text, and if you click on the leaves (the tips) or the nodes (the interior bifurcations) of the tree, you can color code them. IF you wanted to save a picture of your tree, you could click on the “Export” tab, choose the PDF format, and click “Export.” It should pop up in a new tab.


## Post-lab assignment

For your post-lab assignment, answer the critical thinking question below.

**Critical thinking question:**

As with last week, develop a simple question about your project dataset that you can answer using BLAST, making trees, or both.

For example, your question could be something like: 'I found a cytochrome c gene in my dataset. Where does it fit on a tree of other cytochrome C genes from KEGG?' (Make a tree)

Or, "How many ORFs in my dataset have a match to a viral gene? How does that compare to a dataset from a different depth from the same part of the ocean?" (Use BLAST)

*For this week's postlab assignment, describe:*

-What question did you ask?

-How did you go about answering it? (Write this like you would a mini Materials and Methods section: concisely include all the important details needed to repeat what you did, like which databases you searched, which software packages you used, and which important flags you used in your commands. You should include enough information for an intelligent researcher to be able to replicate your results.)

-What were your results? Describe them. Show a table or a figure if appropriate.

-What do you think this tells you about your project dataset? (Think of this as a mini Discussion section: I'm looking for evidence that you thought about your results and how they connect more broadly to some ecological or evolutionary pattern in your dataset.)

Submit via Moodle by lab time next week. **I prefer to grade these blind, so please put your student ID number, rather than your name, on your assignment. (This applies to all future assignments as well.)**
