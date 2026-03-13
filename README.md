### Simulations of the evolution of ecDNA-amplified gene profiles
Simulation code for growing ecDNA-driven populations, tracking oncogene and passenger gene content over time during resistance challenge experiments. Code simulates serial passaging _in vitro_ cell cultures, under successively stronger drug concentration. Series of drug concentrations are modelled by a vector of parameter pairs _m_ and _n_. The values of vector can be changed by editing the nm_pairs vector in the simulation code. Tested using Apple clang version 17.0.0 (clang-1700.6.4.2). Compile with the command

> g++ -o ecDNA_gene_profile_optimisation ecDNA_gene_profile_optimisation.cpp.cpp -std=c++20

and run the simulations using the command

> ./ecDNA_gene_profile_optimisation 
