% Driver script for running map

% Input:
% Sequences (or set of sequences)
% Heterochromatin profile for this set of sequences

noSequences = length(Seq_DNA);
profiles = cell(1,noSequences);

% For each sequence
% Calculate map

for i=1:noSequences
    seq_to_calculate = Seq_DNA(i).sequence;
    seq_to_calculate(Seq_DNA(i).centre_point) = -seq_to_calculate(Seq_DNA(i).centre_point); % flip centre point of sequence
  
    % Uncomment the next four lines if methylation and/or occupancy is required in the simulation
   
    % methylation file ('chrN-start-end_WT.txt'): at base pair accuracy
    % chromosome   start   end   methylation (averaged 0-1)
    %
    % occupancy file ('chrN-start-end_WT.txt'): at whatever resolution (assumped to be 10bp windows)
    % chromosome  start_region  end_region .. ...   chromosome start end  occ1 occ2
    % with the end of the window in column 6 and averaged computed occupancy with the window in cols 7 and 8 for two replicates.

    % methylation = get_methylation(Seq_DNA(i).header, 'WT','../TET_knockout/Methylation/');
    % occupancy = get_occupancy(Seq_DNA(i).header, 'WT',  '../TET_knockout/Nuc_occupancy/',nuc_density,6,7,8);
    
    % methylation = methylation(1:length(seq_to_calculate),2);
    % occupancy = occupancy(1:length(seq_to_calculate),2);
    
    % Assume no methylation or occupancy data

    methylation = [];
    occupancy = [];
    
    file_out = ['WT-' file '-'  Seq_DNA(i).header '.txt'];

    file_out = strrep(file_out, ':', '.');

    disp(file_out);

    % Compute map with any combination of methylation and occupancy (so a 2x2 design if meth & occ given)
    
    profiles{i} = CalculateMapAll(seq_to_calculate,file_out,methylation,occupancy);
end

% write output to file

zz = strrep(Seq_DNA(i).header,':','.');

save(['WT-all-' zz '.mat'],'profiles');

