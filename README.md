ClusterFinder
=============

Predicting biosynthetic gene clusters in genomes
Authors: Peter Cimermancic & Michael Fischbach


Requirements:
 - python (2.X)
 - numpy

Instructions:
 - an example of input file: example_input.txt:
   COLUMN DESCRIPTION:
	1. GeneID
	2. Sequencing status
	3. Organism name
	4. Scaffold OID
	5. Organism OID
	6. Locus Tag
	7. Gene Start
	8. Gene End
	9. Strand
	10. Pfam Template Start
	11. Pfam Template End
	12. Pfam Start
	13. Pfam End
	14. PfamID
	15. Pfam E-score
	16. Enzyme ID
 - if your input file format differs from the one above, please modify the file
   or lines 51-55 of the ClusterFinder.py script

 - an example of running ClusterFinder is shown in ClusterFinder.py script
   DESCRIPTION:
	1. modify paths (if not running from ClusterFinder directory - lines 7, 19 & 20)
	2. name the organism and the input file - lines 15 & 16
	3. run: python ClusterFinder.py
 - testing
   run: python ClusterFinder.py
   without making any changes to the files

 - OUTPUT1 [organims_name.out]: same as input + a column with probability values	
 - OUTPUT2 [organism_name.clusters.out]: same as OUTPUT1, but only for the domains
					from gene clusters that have passed the
					filtering steps.
