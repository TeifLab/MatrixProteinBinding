function [x] = get_occupancy (FASTA_header,file,dir,density,col1,col2,col3)

if isempty(dir)
    dir='.';
end

if(nargin<6)
    col1=14;
    col2=18;
    col3=19;
end

temp = strrep(FASTA_header,':',' '); % to get around string processing below breaking and reading whole

C = textscan(temp,'chr%s %d-%d');

chromosome = char(C{1}{1});

start = C{2};
finish = C{3};

x = zeros(finish-start+1,2);

x(:,1) = start:finish;

y = readtable(strcat(dir,'/',file,'_chr',chromosome,'.',num2str(start),...
    '.',num2str(finish),'.bed'),'ReadVariableNames',false,'Filetype','text');

scaling_1 = 1;
scaling_2 = 1;

if (~isempty(density)) % if nucleosome scalings are provided, otherwise don't scale
    scaling_1 = density(:,{[file '_rep1_scaling']});
    scaling_2 = density(:,{[file '_rep2_scaling']});
    
    try
        row = eval(chromosome); % numeric chromosome name 
    catch
        switch chromosome
            case 'M'
                row = 20;
            case 'X'
                row = 21;
            case 'Y'
                row = 22;
            otherwise
                row = 23;
        end
    end
    
    scaling_1 = table2array(scaling_1(row,1));
    scaling_2 = table2array(scaling_2(row,1));
end

%scaling=0.5*(scaling_1+scaling_2);


try
    y1 = [y(:,col1) y(:,col2)];
    
    y1 = table2array(y1);
    
    y1 = sortrows(y1);
	
	y2 = [y(:,col1) y(:,col3)];

    y2 = table2array(y2);
    
    y2 = sortrows(y2);
    
    for i=start:finish % fill in occupancy using data, note that quoted position
        % is at end of 10bp window for occupancy calculation
        row = i-start+1;
        temp = find(y1(:,1)==(ceil(x(row,1)/10)*10));
        occ1=0;
        if (~isempty(temp))
            occ1=y1(temp(1),2);
        end
        
		temp = find(y2(:,1)==(ceil(x(row,1)/10)*10));
        occ2=0;
        if (~isempty(temp))
            occ2=y2(temp(1),2);
        end
        
        x(row,2)=0.5*(occ1/scaling_1+occ2/scaling_2);
        
    end
    
    
catch
    warning('No occupancy data for this region');
end

end
