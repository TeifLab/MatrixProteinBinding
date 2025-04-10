% Main calculating loop
%
% Calculates binding map (only) based on input sequence
% for CTCF profile and heterochromatin markers

function [profile] = CalculateMap(varargin)

global Lpolymer fNumberOfLigands rank NumberOfPoints eNumberOfChromatinStates R
global X_spacer L_spacer R_spacer
global Seq Seq1 % polymer sequence:  integer, dimension (1:10000)
global m mStart % ligand lengthes: integer, dimension (1:10)
global mm % (the same as m(g) but real, not integer: double precision, dimension (1:10)
global nMaxGap % max interaction lengthes : integer, dimension (1:10)
global MAXnMaxGap MINnMaxGap % max(nMaxGap), min(nMaxGap)
global w %cooperativity parameters : double precision w(0:10000,0:10,0:10, 0:5)
global KKK % double precision(10000,10,147,5)
global Unwrap % double precision Unwrap(0:200,1:10)
global cPolymer c0 % c0(10)
global C3D C2D Cloop Alpha Beta
global productKKK
global s sigma % s(1:10000,1:5), sigma(1:5,1:5)
global iGap % states corresponding to g-j-g gaps integer, dimension (1:10000,0:10)
global iLeftFreeEnd iRightFreeEnd % polymer end matrix state numbers
global ManualC01Range ManualC02Range ManualC03Range ManualC04Range ManualC05Range
global ManualBi CalculateMap CalculateCurves Calculate3DMap CalculateXspacer
global noGaps noLongloops
global lig1mod lig2mod lig3mod lig4mod lig5mod lig6mod

switch nargin
    case 4
        occupancy = varargin{4};
        methylation = varargin{3};
        DNA_seq = varargin{1};
        fileout = varargin{2};
    case 3
        occupancy = [];
        methylation = varargin{3};
        DNA_seq = varargin{1};
        fileout = varargin{2};
    case 2
        occupancy = [];
        methylation = [];
        DNA_seq = varargin{1};
        fileout = varargin{2};
    case 1
        occupancy = [];
        methylation = [];
        fileout = '';
        DNA_seq = varargin{1};
        
    otherwise
        error('Not enough input arguments');
        
end

% Replacements for dialog box stuff:
%global bDefaults
%global Cmin Cmax
%global InputFile OutputFile FastaFile SequenceFile

% I/O text streams
%global fid7 fid8 % streams 7,8 from Fortran codes

SetDefaults; % set default values of parameters (if necessary)

% formats

format1col = '%4s';

format_I1F = '%4d';


nMaxGap = nMaxGap*ones(fNumberOfLigands,1);

if(nMaxGap(1)<=1)
    noGaps=true;
else
    noGaps=false;
end


fid8 = -1;

if ~isempty(fileout)
    fid8=fopen(fileout,'w');
end


meth_flag=false;
occ_flag=false;
state_flag=false;

if(~isempty(methylation))
    meth_flag=true;
end
if(~isempty(occupancy))
    occ_flag=true;
end


%binding map calculation:


% do no meth/occupancy

ParametersInitProteinBinding(DNA_seq,[],[],[]); % fill in K parameters, s, sigma
ConstantsInitProteinBinding; % fill in remaining Unwrap values &c

disp('Without any corrections');

[cMap,tetaMap]=MapOfBindingCalc;

% Now variables are defined, fix headers and format

if(eNumberOfChromatinStates>1)
    state_flag=true;
end


zTotalOutputs=meth_flag+occ_flag +... %data columns
    +fNumberOfLigands+eNumberOfChromatinStates*state_flag + ... % for unaffected 
    +(fNumberOfLigands+eNumberOfChromatinStates*state_flag)*(meth_flag+occ_flag+meth_flag*occ_flag); % if meth and/or occ provided

header=cell(1+zTotalOutputs,1);

header{1}='Site';
j=1;
if(meth_flag)
    header{j+1}='Methylation state';
    j=j+1;
end
if(occ_flag)
    header{j+1}='Av. nuc. occupancy';
    j=j+1;
end
for i=1:fNumberOfLigands
    header{j+1}=strcat('Ci',num2str(i));
    j=j+1;
end
if(state_flag)
    for i=1:eNumberOfChromatinStates
        header{j+1}=strcat('teta',num2str(i));
        j=j+1;
    end
end

if(meth_flag)
    for i=1:fNumberOfLigands
        header{j+1}=strcat('Ci',num2str(i),' (inc meth.)');
        j=j+1;
    end
    if(state_flag)
        for i=1:eNumberOfChromatinStates
            header{j+1}=strcat('teta',num2str(i),' (inc meth.)');
            j=j+1;
        end
    end
end

if(occ_flag)
    for i=1:fNumberOfLigands
        header{j+1}=strcat('Ci',num2str(i),' (inc occ.)');
        j=j+1;
    end
    if(state_flag)
        for i=1:eNumberOfChromatinStates
            header{j+1}=strcat('teta',num2str(i),' (inc occ.)');
            j=j+1;
        end
    end
end

if(meth_flag & occ_flag)
    for i=1:fNumberOfLigands
        header{j+1}=strcat('Ci',num2str(i),' (inc meth.+occ.)');
        j=j+1;
    end
    if(state_flag)
        for i=1:eNumberOfChromatinStates
            header{j+1}=strcat('teta',num2str(i),' (inc meth.+occ.)');
            j=j+1;
        end
    end
end
   


if(meth_flag)
    ParametersInitProteinBinding(DNA_seq,methylation,[],[]); % fill in K parameters, s, sigma
    ConstantsInitProteinBinding; % fill in remaining Unwrap values &c

    disp('With correction for methylation');
    
    [cMap_m,tetaMap_m]=MapOfBindingCalc;
end

if(occ_flag)
    ParametersInitProteinBinding(DNA_seq,[],occupancy,[]); % fill in K parameters, s, sigma
    ConstantsInitProteinBinding; % fill in remaining Unwrap values &c

    disp('With correction for nucleosome occupancy');

    [cMap_o,tetaMap_o]=MapOfBindingCalc;
end

if(occ_flag & meth_flag)
    ParametersInitProteinBinding(DNA_seq,methylation,occupancy,[]); % fill in K parameters, s, sigma
    ConstantsInitProteinBinding; % fill in remaining Unwrap values &c

    disp('With correction for both methylation and nucleosome occupancy');

    [cMap_mo,tetaMap_mo]=MapOfBindingCalc;
end

profile = [methylation'; occupancy'; cMap(2:fNumberOfLigands+1,:)];

if(state_flag)
    profile = [profile; tetaMap(1:eNumberOfChromatinStates,:)];
end

if(meth_flag)
    profile = [profile; cMap_m(2:fNumberOfLigands+1,:)];
    if(state_flag)
        profile = [profile; tetaMap_m(1:eNumberOfChromatinStates,:)];
    end
end

if(occ_flag)
    profile = [profile; cMap_o(2:fNumberOfLigands+1,:)];
    if(state_flag)
        profile = [profile; tetaMap_o(1:eNumberOfChromatinStates,:)];
    end
end

if(occ_flag & meth_flag)
    profile = [profile; cMap_mo(2:fNumberOfLigands+1,:)];
    if(state_flag)
        profile = [profile; tetaMap_mo(1:eNumberOfChromatinStates,:)];
    end
end

if (and(~isempty(fileout),fid8>0))

    format = [format1col repmat(' %22s',1,zTotalOutputs) '\n'];
    fprintf(fid8, format, header{:});
end

formats={format_I1F};


for TestSiteNumber=1:Lpolymer
    
    if (and(~isempty(fileout),fid8>0))
        format = [format_I1F repmat('%23.10E',1,zTotalOutputs) '\n'];
        fprintf(fid8,format, TestSiteNumber, profile(:,TestSiteNumber)');
    end
    
end

if (fid8>0)
    fclose(fid8);
end


end

