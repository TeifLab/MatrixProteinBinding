# MatrixProteinBinding
Uses the transfer matrix method to compute protein binding at single base pair resolution including competition from nucleosomes and DNA methylation if provided. The assembly of the trabnsfer matrices is based on [https://github.com/TeifLab/ChromHL](ChromHL), which originates from the Fortran code by [https://github.com/epigenereg](epigenereg) translated to MATLAB by [https://www.github.com/geejaytee](geejayee), and TFaffinity code by [https://www.github.com/geejaytee](geejayee) translating the R code from the TRAP sopftware from Martin Vingron's lab (https://trap.molgen.mpg.de/cgi-bin/trap_form.cgi).

## License
This software codes belongs to the Laboratory of Gene Regulation at the University of Essex ([generegulation.org](generegulation.org)) and is availabe under CC-BY-ND-NC licence (non-commercial, no derivateves, full attribution). Please provide references to [https://github.com/TeifLab/ChromHL](ChromHL) and the following publication:

Thorn G.J., Clarkson C.T., Rademacher A., Mamayusupova H., Schotta G., Rippe K. and Teif V.B. (2022) DNA sequence-dependent formation of heterochromatin nanodomains. Nature Communications 13, 1861

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

## Subdirectories
- ```MATLAB/``` contains the  MATLAB version by [https://www.github.com/geejaytee](geejayee).
