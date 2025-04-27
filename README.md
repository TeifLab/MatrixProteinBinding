# MatrixProteinBinding
Uses the transfer matrix method to compute protein binding at single base pair resolution including competition from nucleosomes and DNA methylation if provided.

## Inputs required
- Code requires sequences in FASTA format, one per file, with the header:
````
>chrN:start-end
<sequence>
````
- If methylation and nucleosome occupancy provided, they need to be in the format as specified in the Driver.m script, and the competition functions in ParametersInitProteinBinding.m changed to what functional dependence the binding on the methylation or occupancy is modelled.
- If more than one protein to be modelled, then the second and subsequent proteins need to be in the same format as the one given (CTCF), and the coooperativity matrix w needs to change,


## Outputs
A fixed-width format file given the binding probability for each protein per nucleotide along the sequence: if DNA methylation and/or occupancy are provided then binding probabilities for all combinations of methylation and nucleosome occupancy are calculated too.

## Notes
Based on a FORTRAN code by  [https://github.com/epigenereg](Vladimir Teif) , using sections of ChromHL, for the assembly of the transfer matrix, and TFaffinity, for the calculation of the TRAP-derived binding affinity used as input.

## Subdirectories
- ```MATLAB/``` contains the original MATLAB version by [https://www.github.com/geejaytee](geejayee).
