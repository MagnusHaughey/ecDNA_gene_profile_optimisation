### Simulations of the evolution of ecDNA-amplified gene profiles
Simulation code for growing ecDNA-driven populations, tracking oncogene and passenger gene content over time during resistance challenge experiments. Code simulates serial passaging _in vitro_ cell cultures, under successively stronger drug concentration. Series of drug concentrations are modelled by a vector of parameter pairs _m_ and _n_. The values of vector can be changed by editing the nm_pairs vector in the simulation code. Tested using Apple clang version 17.0.0 (clang-1700.6.4.2). Compile with the command

```
g++ -o ecDNA_gene_profile_optimisation ecDNA_gene_profile_optimisation.cpp.cpp -std=c++20
```

and run the simulations using the command

```
./ecDNA_gene_profile_optimisation --verbose  -N [N] -k [k] -s [s] -l [l] -p [p] -q [q] -x [x]
```

where\
&nbsp; --verbose &emsp;&emsp; verbose flag (optional)\
&nbsp; -N &emsp;&emsp; maximum population size (in number of cells)\
&nbsp; -k &emsp;&emsp; ecDNA copy number in initial cell\
&nbsp; -s &emsp;&emsp; selection strength. Scalar multiplier for the ecDNA dependent birth rate function, with s=0 giving rise to neutral growth and s>0 giving rise to positive ecDNA copy number dependent selection\
&nbsp; -l &emsp;&emsp; ecDNA gene burden penalty\
&nbsp; -p &emsp;&emsp; ecDNA fusion rate\
&nbsp; -q &emsp;&emsp; ecDNA splitting rate\
&nbsp; -x &emsp;&emsp; random seed

Flags should be specified before numerical arguments.

Data are written a directory named ./RESULTS/Nmax=`N`_resampleSize=20000_k=`k`_s=`s`_l=`l`_p=`p`_q=`q`_nsequence=`n1_n2_n3...`_msequence=`m1_m2_m3...`_seed=`x`/n=`n`_m=`m`_tissue.csv
