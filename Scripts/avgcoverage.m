% Calculate average coverage: count the number of reads that align to a
% locus for all samples, then divide by the number of samples (79). With
% this coverage value for each locus, find the average across all 3015
% loci.

load('t.mat');
C = zeros(size(t.sequences,1),1);

for n = 1:size(t.sequences,1)
    cd '/Users/jonathanchan/Documents/MATLAB/epigenetics'
    template = t.sequences{n};
    locus = t.loci{n};

    if size(find(template == 'G')) == [1,0]
       fprintf(2,'\nNo cpg sites at locus %s \n', locus);
    else
        cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/bam';
        myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/bam';
        bamFiles = dir(fullfile(myDir, '*.bam'));

        % find how many reads are in locus n (across all samples)
        locus_n_count = 0;
        for m = 1:size(bamFiles,1)
            fullName = fullfile(bamFiles(m).name);
            command = sprintf('samtools view %s %s | cut -f10',fullName,locus);
            [status,cmdout]=unix(command);
            reads = regexp(cmdout, '\n', 'split');
            locus_n_count = locus_n_count + size(reads,2);
        end

        C(n) = locus_n_count/79; % mean(C) = 326.9116
    end
end