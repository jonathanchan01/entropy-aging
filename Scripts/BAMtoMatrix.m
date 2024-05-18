% Given a single locus and the corresponding Hg38 template, return CHALM, 
% average methylation, and entropy values

function metrics = BAMtoMatrix(template, locus, reverse)

% set directory and pull bam files
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data/addtl_bam';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data/addtl_bam';
addpath '/Users/jonathanchan/Documents/MATLAB/epigenetics/Scripts'
bamFiles = dir(fullfile(myDir, '*.bam'));

% hold all the relevant data for each bam file
m = struct('file',{},'cpgSites', {}, 'alignment', {}, 'matrix',{},'avgMethylation',{});
metrics = struct('ch',{},'aa',{},'en',{},'age',{});
fprintf(1, '\nReading from %s\n', locus)

% loop through all bam files
for n = 1:size(bamFiles,1)
    
    fullName = fullfile(bamFiles(n).name);
    command = sprintf('samtools view %s %s | cut -f10',fullName,locus);
    [status,cmdout]=unix(command);
    reads = regexp(cmdout, '\n', 'split');
    m(n).file = bamFiles(n).name;
    metrics(n).age = sampleAges(bamFiles(n).name);
    
    % output progress
    switch(n)
        case ceil(size(bamFiles,1)/4)
            fprintf(2,'25%% done\n');
        case ceil(size(bamFiles,1)/2)
            fprintf(2,'50%% done\n');
        case ceil(size(bamFiles,1)/4*3)
            fprintf(2,'75%% done\n');
    end
    
    % only keep data for bam files that have data at that locus
    if size(reads,2) <= 2
        fprintf(1,'Not enough data at %s from sample %s\n',locus,m(n).file);
    else
        fasta = makeFasta(reads,template);
        aligned = alignSequences(fasta);
        cpgSites = find(aligned(1, :) == 'G');        
        result = makeBinMatrix(cpgSites,aligned,reverse);
        
        % store data
        m(n).cpgSites = cpgSites;
        m(n).alignment = aligned;
        m(n).matrix = result;
        m(n).avgMethylation = avgMeth(result);
    end

end

% calculate metrics
CHALMresult = CHALM(m);
for c = 1:size(CHALMresult,1)
    metrics(c).ch = CHALMresult(c);
end

AAresult = AA(m);
for a = 1:size(AAresult,1)
    metrics(a).aa = AAresult(a);
end

[etpyM, sampleIDs_counts, m_size] = makeEtpyMatrix(m);
ENresult = calcEtpy(etpyM,sampleIDs_counts, m_size);

for e = 1:size(ENresult)
    metrics(e).en = ENresult(e);
end

fprintf(1, 'Finished reading %s\n', locus);
end