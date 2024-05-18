% make binary matrix of methylation data, where each row represents a read
% and each column represents a CpG site. 1 = methylated, 0 = unmethylated.

function binMatrix = makeBinMatrix(cpgSites,aligned,reverse)

binMatrix = zeros(size(aligned,1), nnz(cpgSites));

% G = 1, A = 0 (data from reverse strand)
% everything else = NaN
for c = 1:nnz(cpgSites)
    currentSite = cpgSites(c);
    for r = 1:size(aligned, 1)
        if reverse
            if aligned(r, currentSite) == 'G'
                binMatrix(r,c) = 1;
            elseif aligned(r, currentSite) ~= 'A'
                binMatrix(r,c) = NaN;
            end
        else
            if aligned(r, currentSite) == 'C'
                binMatrix(r,c) = 1;
            elseif aligned(r, currentSite) ~= 'T'
                binMatrix(r,c) = NaN;
            end
        end
    end
end

% delete the row for the reference sequence
binMatrix(1,:) = [];
