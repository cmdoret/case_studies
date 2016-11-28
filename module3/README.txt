
#########
Directories
#########

PDZLigands/: fasta files of all phage peptide for each domain. The file title indicates the 'sourceProtein'-'whichPDZdomain'
PDZLigands_LOLA: files in the LOLA format for all phage ligands for each domain. The project.txt files should be loaded in LOLA to display all logos

lola-1.1/ software for visualizing the logos. For Mac, double click lola-1.1-beta.jar file. For Windows, double click lola.bat. On Linux run lola.sh. 

#########
Files
#########

Alignment files:
******The first entry shows the position of the binding site ('B').*********
PDZ_SMART_CLUSTAL_sub.fa: Alignment of PDZ sequences retrieved from SMART database. The alignment was done on all human and worm PDZ sequences with CLUSTAL_OMEGA
PDZ_SMART_MUSCLE_sub.fa: Alignment of PDZ sequences retrieved from SMART database. The alignment was done on all human and worm PDZ sequences with MUSCLE
PDZ_phage_MUSCLE.fa: Alignment of PDZ sequences used in the clones for the phage experiment. The alignment was done on human and worm PDZ sequences considered here with MUSCLE

Codon bias:
phageLibraryNNKTheoreticalCodonBias.txt: file containing the theoretical codon bias (expected frequency of each amino acid in a random phage library).

Paper:
Tonikian.pdf: Original paper (Supplemetary Figures available online at http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060239).

Instructions:
Instructions.pdf: detailed instructions on how to proceed to reproduce the Figures.

PDZ specificity classes:
PDZclass.txt: File containing different classification schemes. Second column: The 2 main classes. Third column: 7 classes (manually defined). Fourth column: the classes defined in Tonikian et al.

analyze.R: Script to load the data. You should fill the holes in this code to perform the analysis.
