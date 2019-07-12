tic;

nuc_density = readtable('Nucleosome densities.txt');

system(['mkdir -p Output_files/; mkdir -p Text_output/']);

% Adjust the following line for the FASTA files which should be one per region and be of the form
% >chrN:start-end
% <sequence> - either upper or lower case

fasta_dir = 'common_FASTA/'


try
    Seq_DNA=Read_FASTA_all([ fasta_dir '/*.fa'])
    if ~isempty(Seq_DNA)
        Driver;
        system(['mv *chr*.txt Text_output/']);
    end
end
catch message
    display(['ERROR in file: ' message.stack.file])
    display(['ERROR: ' getReport(message)])
end

toc;

exit

