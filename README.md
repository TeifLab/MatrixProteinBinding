# MatrixProteinBinding (part of ChromHL)
Uses the transfer matrix method to compute protein (e.g. transcription factor, TF) binding to DNA at single base pair resolution considering the competition of TF binding with nucleosomes and the effect of DNA methylation (optional). This is part of the ChromHL software ([https://github.com/TeifLab/ChromHL](ChromHL)) and the same license restrictions apply (CC-BY-NC-ND). The assembly of the transfer matrices is based on the ChromHL code ([https://github.com/TeifLab/ChromHL](ChromHL)), which originates from the Fortran code by [https://github.com/epigenereg](epigenereg) translated to MATLAB by [https://www.github.com/geejaytee](geejayee). The assignment of transcription factor affinity is based on the TFaffinity ([https://github.com/TeifLab/TFaffinity](TFaffinity)) code by [https://www.github.com/geejaytee](geejayee) which is using the TRAP algorithm (Martin Vingron's lab, https://trap.molgen.mpg.de/cgi-bin/trap_form.cgi).

## License
This software belongs to the Laboratory of Gene Regulation at the University of Essex ([generegulation.org](generegulation.org)) and is availabe to external users under CC-BY-NC-ND licence (non-commercial use only, no derivateves, full attribution required, including the name of the software ChromHL). Please provide references to [https://github.com/TeifLab/ChromHL](ChromHL) and the following publication:

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
