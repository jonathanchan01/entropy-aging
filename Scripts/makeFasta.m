% create fasta file

function fasta = makeFasta(reads,template)

% first row is the template sequence:
% bisulfite conversion: G -> A except when adjacent to C
temp=fopen('fasta.txt','w');
fprintf(temp,'>alignment\n');
fprintf(temp,template);
fprintf(temp,'\n');
for k=1:size(reads,2)-1
    fprintf(temp,'>read %d\n',k);
    fprintf(temp, '%s\n', reads{k});
end

fasta = fastaread('fasta.txt');
fclose(temp);

end