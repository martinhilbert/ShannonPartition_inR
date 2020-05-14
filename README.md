# ShannonPartition_inR
Shannon information partition from empirical data and from distribution
The script calculates the entropy and mutual information for all possible relations from several variables. It calcualtes all "areas" and intersections of an information-theoretic Venn diagram (called information diagram)
https://en.wikipedia.org/wiki/Multivariate_mutual_information

The script either reads in empirical data (such as a time sequence, from which it derives the underlying distribution) or works directly with a predefined probability distribution. 
It was adopted from an earlier Python implementation by Ryan James: https://github.com/Autoplectic/dit 

The final output prints all partitions + the user is able to print out different subsets of sums of the partition, i.e. for the conditional entropies. For example: H(A | B) would be the sum of all atoms that quantify the entropy of A and are conditioned on B. For example, for three variables: H(A|B) = H(A|B;C) + I(A;C|B)

The scrip called "code.R" is the main script. Open it in R. The other scripts (and data) simply need to be placed within the defined working directory
