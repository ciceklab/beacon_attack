# beacon_attack
The repository contains three files that are three main functions (QI_Method, Lambda_QI and GI_Method).

QI_Method.m contains the code to query a simulated beacon with the SNP information available. The received answer is then used in
Lambda_QI.m to determine the Lambda values for that individual and the QI attack. The GI_Method.m computes the responses the beacon
gives to the GI attack queries. The response set can then be used to calculate the Lambda values for the individuals as shown in
the paper.

The data we used is available on HapMap (ftp://ftp.ncbi.nlm.nih.gov/hapmap).

In order to query the real beacons we used the beacon-network.org interface and data from the personal genomes project (https://my.pgp-hms.org/public_genetic_data).
