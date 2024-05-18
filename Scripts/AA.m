% For each binary matrix, return the average methylation. 

% For each sample, create an n x 1 vector for each read. 
% Store the average methylation of each vector of the sample into vector v 
% of averages. The result for each sample is the average of v.

function AAresult = AA(m)

AAresult = [];

%loop through m
for nn = 1:size(m,2)
    v = [];
    
    % at each row...
    for r = 1:size(m(nn).matrix,1)
        count1 = 0;
        count = 0;
        
        % calculate the avg methylation
        for c = 1:size(m(nn).matrix,2)
            if m(nn).matrix(r,c) == 1
               count1 = count1+1;
               count = count+1;
            end
            if m(nn).matrix(r,c) == 0
               count = count+1; 
            end
        end
        
        a = count1/count;
        
        if isnan(a)
           a = 0; 
        end
        
        v = [v; a];
    end
    
    AAresult = [AAresult; mean(v)];
    
end