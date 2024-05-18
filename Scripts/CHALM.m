% For each binary matrix, return the value of the CHALM metric. 

% Create a binary n x 1 vector v for each sample, where n = number of reads.
% If any read i in the row is methylated, v_i = 1. If not, v_i = 0. The 
% result for each sample is the average of v.

function CHALMresult = CHALM(m)

CHALMresult = [];

%loop through m
for nn = 1:size(m,2)
    v = [];
    
    % loop through matrix nn
    for r = 1:size(m(nn).matrix,1)
        for c = 1:size(m(nn).matrix,2)
            % create binary vector; 1 if any read in the row is methylated
            if m(nn).matrix(r,c) == 1
                v = [v; 1];
                break;
            end
            if c == size(m(nn).matrix,2)
                v = [v; 0];
            end
        end
    end
    
    CHALMresult = [CHALMresult; mean(v)];
    
end