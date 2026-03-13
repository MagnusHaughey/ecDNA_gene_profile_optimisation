### Simulations of the evolution of ecDNA-amplified gene profiles
Simulation code for growing ecDNA-driven populations, tracking oncogene and passenger gene content over time during resistance challenge experiments. Code simulates serial passaging _in vitro_ cell cultures, under successively stronger drug concentration. Series of drug concentrations are modelled by a vector of parameter pairs `m` and `n`. The values of vector can be changed by editing the `nm_pairs` vector in the simulation code. Tested using Apple clang version 17.0.0 (clang-1700.6.4.2) and Python 3.9.6. Compile with the command

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

Data are written a directory named ./RESULTS/Nmax=`N`_resampleSize=20000_k=`k`_s=`s`_l=`l`_p=`p`_q=`q`_nsequence=`n1_n2_n3...`_msequence=`m1_m2_m3...`_seed=`x`/n=`ni`_m=`mi`_tissue.csv

where `n1_n2_n3...` and `m1_m2_m3...` are the sequences of n and m values specified in the `nm_pairs` vector, and `ni` and `mi` are the `n` and `m` values of the serial passage to which the data relates.


Compute and plot gene colocalization and copy number data for the resulting files by executing

```
python3 make_plots_increasingDrugConcentration.py
```
Manually edit the `all_n_sequences` list in `make_plots_increasingDrugConcentration.py` to match the desired nsequence of the raw data in ./RESULTS/. The code will iterate over all parameter combination corresponding to this sequence and plot colocalization and gene copy number data for each combination. Data is plotted in a new file in the ./PLOTS/ directory
