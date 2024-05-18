% align sequences

function aligned = alignSequences(fasta)

% only use the bases that aligned well
for j = 1:size(fasta,1)
    [Score, Alignment, Start] = swalign(fasta(1,:),fasta(j,:),'Alphabet','NT');
    read = Alignment(3,1:size(Alignment,2));
    final = erase(read,'-');
    fasta(j).Sequence = final;
end

seqs = {fasta.Sequence};

% align
aligned = multialign(seqs, 'GapOpen', 50, 'terminalGapAdjust',true);