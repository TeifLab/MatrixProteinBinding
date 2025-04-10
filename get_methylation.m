function [x] = get_methylation (FASTA_header,file,dir)
 
if nargin==2
    dir='.';
end
 
temp = strrep(FASTA_header,':',' '); % to get around string processing below breaking and reading whole 
 
C = textscan(temp,'chr%s %d-%d');
 
chromosome = char(C{1}{1});
 
start = C{2};
finish = C{3}; 
 
x = zeros(finish-start+1,2);
 
x(:,1) = start:finish;
 
y = readtable(strcat(dir,'/',file,'_chr',chromosome,'_',num2str(start),...
    '_',num2str(finish),'_meth.txt'),'ReadVariableNames',false);
 
y = table2array(y(:,2:end));
 
y = sortrows(y);

try 
  x(ismember(x(:,1),y(:,1)),2)=y(:,4);
catch
   disp('No methylation data for this region');
end
