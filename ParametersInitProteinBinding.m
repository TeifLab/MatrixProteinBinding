function [] = ParametersInitProteinBinding(DNA_seq,methylation,occupancy,func)

% input: DNA_seq (per base)

global Lpolymer fNumberOfLigands  eNumberOfChromatinStates
global m % ligand lengthes: integer, dimension (1:10)
global mm % (the same as m(g) but real, not integer: double precision, dimension (1:10)
global nMaxGap % max interaction lengthes : integer, dimension (1:10)
global w %cooperativity parameters : double precision w(0:10000,0:10,0:10, 0:5)
global KKK % double precision(10000,10,147,5)
global Unwrap % double precision Unwrap(0:200,1:10)
global c0 % c0(10)
global K0 % K0(10) % prefactor for PWM
global s sigma % s(1:10000,1:5), sigma(1:5,1:5)
global noGaps
global CTCF_profile

fNumberOfLigands=1; % CTCF
eNumberOfChromatinStates=1; % no need for other states

if (isempty(func))
    func = @(x) exponential(x,log(15));  %
end 

Lpolymer = length(DNA_seq);

s = zeros(Lpolymer,5);
sigma = zeros(5,5);
Unwrap = zeros(2001,10); % note that first component was zero-based in F90 code
KKK = zeros(Lpolymer,10,20,5);

% First ligand (CTCF)
g=1;
m(g)=19; % CTCF binds to 19 nucleotides

mm(g)=real(m(g));
Unwrap(:,g)=0.; %This line should be always equal to zero
Unwrap(1,g)=1.; %This line should be always equal to one 

% concentration
c0(g)=10^-6;

% build PWM affinity curve
K0(g)=10^9;

CTCF_profile = ReadBindingConstants(DNA_seq,'CTCF_matrix_Orlov.txt')); % build CTCF affinity curve


if(isempty(occupancy))
    occupancy = zeros(Lpolymer,1);
end

if(isempty(methylation))
    methylation = zeros(Lpolymer,1);
end


K_macro = CTCF_profile(1:(Lpolymer-m(g)+1));

% binding constants based on PWM value generated from ReadBindingConstants
% v3 code 2017-10-16: corrections are now per base
%
% KKK(i,...,j,...) = 0 if i<j
%                  = K_macro (i-j+1)^(1/n) * (1-meth(i)) * func(occ(i)/m(g))
%                  = 0 if i>length-m(g)+1
%
% Upshot is that binding is corrected by average occupancy over binding site
% and by product of (1-meth) over binding site, so two @ 50% methylation
% is a reduction of 75% ( (1-0.5)*(1-0.5) = 0.25 = (1-0.75) )
%
%%
%%
%- v4 code
%
%- precorrect K_macro by func(av(occ)) before this loop here:

occ = movmean(occupancy,[0 m(g)-1]); % sliding window average occupancy over binding sites in sequence

occ = occ(1:length(K_macro));

K_macro_occ = K_macro .* func(occ)';

% Methylation correction - no binding *if any methylation>0*

correct_methylation = 1-(methylation>0);

for i=1:m(g)

    indices = (1:length(K_macro)) + (i-1);                     % indices of non-zero KKK values
 
    meth_corr = correct_methylation(indices);                  % sliding-window of methylation_correction

    KKK(indices,g,i,1) = K_macro_occ' .^ (1/m(g)) ...          % pre-occ-corrected CTCF profile
                         .* meth_corr;                         % now corrected by methylation_correction


end


% Other proteins go here: copy the above code, change the multiplier for the affinity (K0(g)), 
% change the PWM matrix 'CTCF_matrix_Orlov.txt' to the PWM of that protein
% change the length of the protein m(g)
% change corrections for methylation and/or nucleosome occupancy


% parameters for chromatin states: s
s(:,:)=1.;

% sigma
sigma(:,:)=1.;

% Cooperativity constants (note indices incremented by 1)
w(:,:,:,:)=1.;
% change this w to get different cooperativities between proteins


nMaxGap(1:fNumberOfLigands)=0;
noGaps=false;

end  % subroutine ParametersInitUnwrap
