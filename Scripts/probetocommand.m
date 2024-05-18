% From the txt file of probes, create .bed files for all probes and extract their loci and sequence, using Hg38 as a reference. Return the loci and the corresponding template sequence.
function t = probetocommand()

%% Create .bed files for each probe
% read in probe file
fileID = fopen('human_probes.txt');
temp = textscan(fileID,'%s %s %d %d');
fclose(fileID);

% get only the lines in temp that have the code 'R', 'C', or 'O'
temp_probes.code = temp{1,1};
temp_probes.chr = temp{1,2};
temp_probes.start = temp{1,3};
temp_probes.end = temp{1,4};


% keep only lines that start with 'O'
index = 1;
probes = struct('code',{},'chr',{},'chr_alt',{},'start',{},'end',{});
for n = 1:size(temp_probes.code,1)
    if temp_probes.code{n,1} == 'O'
        probes(index).code = temp_probes.code{n,1};
        probes(index).chr = temp_probes.chr{n,1};
        probes(index).start = temp_probes.start(n,1);
        probes(index).end = temp_probes.end(n,1);
        index = index+1;
    end
end

templatefile = fopen('templates.txt','w+');

% convert chromosome number into RefSeq name
for n = 1:size(probes,2)
   chr = probes(n).chr;
   switch chr
       case '1'
           probes(n).chr_alt = 'NC_000001.11';
       case '2'
           probes(n).chr_alt = 'NC_000002.12';
       case '3'
           probes(n).chr_alt = 'NC_000003.12';
       case '4'
           probes(n).chr_alt = 'NC_000004.12';
       case '5'
           probes(n).chr_alt = 'NC_000005.10';
       case '6'
           probes(n).chr_alt = 'NC_000006.12';
       case '7'
           probes(n).chr_alt = 'NC_000007.14';
       case '8'
           probes(n).chr_alt = 'NC_000008.11';
       case '9'
           probes(n).chr_alt = 'NC_000009.12';
       case '10'
           probes(n).chr_alt = 'NC_000010.11';
       case '11'
           probes(n).chr_alt = 'NC_000011.10';
       case '12'
           probes(n).chr_alt = 'NC_000012.12';
       case '13'
           probes(n).chr_alt = 'NC_000013.11';
       case '14'
           probes(n).chr_alt = 'NC_000014.9';
       case '15'
           probes(n).chr_alt = 'NC_000015.10';
       case '16'
           probes(n).chr_alt = 'NC_000016.10';
       case '17'
           probes(n).chr_alt = 'NC_000017.11';
       case '18'
           probes(n).chr_alt = 'NC_000018.10';
       case '19'
           probes(n).chr_alt = 'NC_000019.10';
       case '20'
           probes(n).chr_alt = 'NC_000020.11';
       case '21'
           probes(n).chr_alt = 'NC_000021.9';
       case '22'
           probes(n).chr_alt = 'NC_000022.11';
       case 'X'
           probes(n).chr_alt = 'NC_000023.11';
       case 'Y'
           probes(n).chr_alt = 'NC_000024.10';
   end
   
   % create .bed files
   bedtool_content = sprintf('%s\t%s\t%s', probes(n).chr_alt, num2str(probes(n).start), num2str(probes(n).end));
   filename = [num2str(n),'_probe.bed'];
   fid = fopen(filename,'w');
   fprintf(fid,bedtool_content);
   
   % get template
   unixcommand = sprintf('bedtools getfasta -fi GCF_000001405.40_GRCh38.p14_genomic.fna -bed %s',filename);
   [status,cmdout]=unix(unixcommand);
   
   % print to file
   fprintf(templatefile, cmdout);
   
   fclose(fid);
end

% Create a vector of template sequences
templatefile = fopen('templates.txt');
temp2 = textscan(templatefile,'%s');
templates = temp2{1,1};
s = size(templates,1);
if mod(size(templates,1),2) == 1
    s = (s+1)/2;
else
    s = s/2;
end

t.sequences = cell(s,1);

for m = 1:s
    % Make all caps
    % t.sequences{m} = upper(templates{2*m});

    % Convert for bisulfite-seq data: if a G does not follow a C, turn it to A
    read = upper(templates{2*m});
    for index = 1:size(read,2)
       if read(index) == 'G'
           % G at the start
           if index == 1
              read(index) = 'A'; 
           end
           
           % G follows something other than C
           if index > 1 && read(index-1) ~= 'C'
               read(index) = 'A';
           end
       end
       t.sequences{m} = read;
    end
    
end


% Create a corresponding vector of loci
t.loci = cell(s,1);

for m = 1:s
    t.loci{m} = strcat('chr',probes(m).chr,':',num2str(probes(m).start),'-',num2str(probes(m).end));
end

fclose(templatefile);

delete *.bed;
end