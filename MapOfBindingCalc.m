%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subroutine MapOfBindingCalc calculates the map ofbinding
% for a given c0(g), g=1,2,3
% ...........dStatda to be added later..............
%
% Stat - statsum
% max - after Trace of the current statsum value reaches max, StatSum
%       is being normilized by max
% NormCount - number of normalizations
% noLeftOverhang=.true. prohibits ligand overhang from the left DNA end
% noRightOverhang=.true. prohibits ligand overhang from the right DNA end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[cMap,tetaMap] = MapOfBindingCalc()

global Lpolymer fNumberOfLigands eNumberOfChromatinStates
global fid7

cMap=zeros([fNumberOfLigands+1  Lpolymer]);
tetaMap=zeros([eNumberOfChromatinStates  Lpolymer]);

TIME=tic;

for TestSiteNumber=1:Lpolymer
    
    fprintf('(site,g,e)={%d,%s,%s}\n',TestSiteNumber,'all','all');
    
    [c,teta]=PointOfMapOfBindingCalc(TestSiteNumber);
    
    cMap(2:end,TestSiteNumber)=c(2:end);
    tetaMap(:,TestSiteNumber)=teta(:);
    
    
end % TestSiteNumber



end

