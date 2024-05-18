% given a single locus and the corresponding Hg38 template, return CHALM, average of
% averages, and entropy values for one sample

function metrics = BAMtoMatrix_customSample(template, locus, file)

% set directory and pull bam files
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/bam';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/bam';
bamFile = dir(fullfile(myDir, file));

% hold all the relevant data for each bam file
m = struct('file',{file},'cpgSites', {}, 'alignment', {}, 'matrix',{},'avgMethylation',{});
metrics = struct('ch',{},'aa',{},'en',{},'age',{},'fileName',{});

fprintf(1, '\nReading from %s\n', locus)

% analyze the single bam file
fullName = fullfile(bamFile(1).name);
command = sprintf('samtools view %s %s | cut -f10',fullName,locus);
[status,cmdout]=unix(command);
reads = regexp(cmdout, '\n', 'split');

m(1).file = bamFile(1).name;
metrics(1).fileName = m(1).file;


switch(m(1).file)
    case 'P11111111111111.bam'
        metrics(1).age = 24.4;
    case 'P11111111111112.bam'
        metrics(1).age = 21.8;
    case 'P11111111111113.bam'
        metrics(1).age = 31.59;
    case 'P41019061105623.bam'
        metrics(1).age = 31.59;
    case 'P41019061105624.bam'
        metrics(1).age = 57.4;
    case 'P41019061105630.bam'
        metrics(1).age = 49.3;
    case 'P41019061105632.bam'
        metrics(1).age = 52.6;
    case 'P41019061105633.bam'
        metrics(1).age = 53.1;
    case 'P41019061105634.bam'
        metrics(1).age = 53.2;
    case 'P41019061105638.bam'
        metrics(1).age = 51.1;
    case 'P41019061105642.bam'
        metrics(1).age = 60.3;
    case 'P41019061105644.bam'
        metrics(1).age = 30.6;
    case 'P41019061105645.bam'
        metrics(1).age = 34.5;
    case 'P41019061105647.bam'
        metrics(1).age = 55.5;
    case 'P41019061105649.bam'
        metrics(1).age = 50;
    case 'P41019061105650.bam'
        metrics(1).age = 33.4;
    case 'P41019061105651.bam'
        metrics(1).age = 64.1;
    case 'P41019061105653.bam'
        metrics(1).age = 47.2;
    case 'P41019061105655.bam'
        metrics(1).age = 40.9;
    case 'P41019061105656.bam'
        metrics(1).age = 60.2;
    case 'P41019061105657.bam'
        metrics(1).age = 47.8;
    case 'P41019061105682.bam'
        metrics(1).age = 21;
    case 'P41019061105686.bam'
        metrics(1).age = 18.5;
    case 'P41019061105688.bam'
        metrics(1).age = 15.1;
    case 'P41019061105693.bam'
        metrics(1).age = 48.4;
    case 'P41019061105696.bam'
        metrics(1).age = 50.8;
    case 'P41019061105701.bam'
        metrics(1).age = 84.4;
    case 'P41019061105735.bam'
        metrics(1).age = 50.3;
    case 'P41019061105737.bam'
        metrics(1).age = 53.2;
    case 'P41019081703384.bam'
        metrics(1).age = 41.6;
    case 'P41019081703391.bam'
        metrics(1).age = 55.4;
    case 'P41019081703456.bam'
        metrics(1).age = 57.7;
    case 'P41019081703461.bam'
        metrics(1).age = 17.9;
    case 'P41019081703462.bam'
        metrics(1).age = 36.1;
    case 'P41019081703465.bam'
        metrics(1).age = 78.2;
    case 'P41019081703471.bam'
        metrics(1).age = 41.8;
    case 'P41019081703479.bam'
        metrics(1).age = 78.2;
    case 'P41019081703500.bam'
        metrics(1).age = 52;
    case 'P41019081703512.bam'
        metrics(1).age = 32.2;
    case 'P41019081703517.bam'
        metrics(1).age = 64.2;
    case 'P41019081703523.bam'
        metrics(1).age = 26.3;
    case 'P41019081703526.bam'
        metrics(1).age = 56.2;
    case 'P41019081703527.bam'
        metrics(1).age = 24.2;
    case 'P41019081703539.bam'
        metrics(1).age = 57.3;
    case 'P41019081703556.bam'
        metrics(1).age = 16.8;
    case 'P41019081703557.bam'
        metrics(1).age = 51.1;
    case 'P41019081703571.bam'
        metrics(1).age = 22;
    case 'P41019081703575.bam'
        metrics(1).age = 49.9;
    case 'P41019081703576.bam'
        metrics(1).age = 41.5;
    case 'P41019081706684.bam'
        metrics(1).age = 73.3;
    case 'P41019081706685.bam'
        metrics(1).age = 37.5;
    case 'P41019081706687.bam'
        metrics(1).age = 54;
    case 'P41019081706707.bam'
        metrics(1).age = 67.8;
    case 'P41019081706739.bam'
        metrics(1).age = 27.2;
    case 'P41019081706743.bam'
        metrics(1).age = 41.8;
    case 'P41019081706745.bam'
        metrics(1).age = 47.9;
    case 'P41019081706746.bam'
        metrics(1).age = 57.7;
    case 'P41019081706759.bam'
        metrics(1).age = 71.3;
    case 'P41019081706782.bam'
        metrics(1).age = 55.7;
    case 'P41019081706783.bam'
        metrics(1).age = 7.2;
    case 'P41019081706785.bam'
        metrics(1).age = 56.1;
    case 'P41019081706789.bam'
        metrics(1).age = NaN;
    case 'P41019081706790.bam'
        metrics(1).age = 50;
    case 'P41019081706814.bam'
        metrics(1).age = 40.7;
    case 'P41019081706836.bam'
        metrics(1).age = 46.4;
    case 'P41019091104671.bam'
        metrics(1).age = 41.8;
    case 'P41019091104677.bam'
        metrics(1).age = 40.7;
    case 'P41019091104727.bam'
        metrics(1).age = 33.5;
    case 'P41019091104730.bam'
        metrics(1).age = 75.6;
    case 'P41019091104731.bam'
        metrics(1).age = 35.2;
    case 'P41019091104746.bam'
        metrics(1).age = 35;
    case 'P41019091104751.bam'
        metrics(1).age = 69.2;
    case 'P41019091104753.bam'
        metrics(1).age = 29.2;
    case 'P41019091104756.bam'
        metrics(1).age = 32.49;
    case 'P41019091109805.bam'
        metrics(1).age = 67.4;
    case 'P41019091109813.bam'
        metrics(1).age = 32.79;
    case 'P41019091109817.bam'
        metrics(1).age = 20.3;
    case 'P41019091109831.bam'
        metrics(1).age = 51.6;
    case 'P41019091109837.bam'
        metrics(1).age = 66;
    case 'P41019091109844.bam'
        metrics(1).age = 30.6;
    otherwise
        fprintf('Error: sample %s does not have an age.',m(1).file);
end

% only keep data for bam files that have data at that locus
if size(reads,2) <= 2
    fprintf(1,'Not enough data at %s from sample %s\n',locus,m(1).file);
else
    cd '/Users/jonathanchan/Documents/MATLAB/epigenetics';
    fasta = makeFasta(reads,template);
    aligned = alignSequences(fasta);
    cpgSites = find(aligned(1, :) == 'G');        
    result = makeBinMatrix(cpgSites, aligned);

    % store data
    m(1).cpgSites = cpgSites;
    m(1).alignment = aligned;
    m(1).matrix = result;
    m(1).avgMethylation = avgMeth(result);
    cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/bam';
end

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics';

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

% calcEtpy_focus
age = metrics(1).age;
ENresult = calcEtpy_focus(etpyM,sampleIDs_counts, m_size, age);

for e = 1:size(ENresult)
    metrics(e).en = ENresult(e);
end


fprintf(1, 'Finished reading %s\n', locus);

end