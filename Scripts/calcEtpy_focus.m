% Calculate the Shannon entropy of rows of the matrix.

% Find all possible entropy states. For each sample with at least 50 reads,
% count the number of times each state is observed in the matrix. Normalize
% the resultant vector to get a vector of probabilities. Using this, 
% calculate the informational entropy using the Shannon entropy formula.
% Return a vector of entropies, one value per sample.


function etpy = calcEtpy_focus(etpyM, countSampleIDs, m_size, age)

% if no data
if size(etpyM) == [0,0]
   etpy = NaN(m_size,1);
   return
end

% find how many CPG sites of interest
n = size(etpyM,2)-1;

% there are 2^n possible states, create a matrix of all possibilities
states = [];
if n > 1
    for l = 0:n
        M = [zeros(1,l),ones(1,n-l)];
        n = numel(M);
        k = sum(M);
        c = nchoosek(1:n,k);
        m = size(c,1);
        binary = zeros(m,n);
        binary(sub2ind([m,n],(1:m)'*ones(1,n-l),c)) = 1;
        states = [states; binary];
    end
end

if n == 1
    states = [0; 1];
end

index = 0;
etpy = [];
for sampleID = 1:size(countSampleIDs,1)
    numReads = countSampleIDs(sampleID);
    sample_etpyM = [];
    index = index+1;
    for index = index:numReads+index-1
        % construct matrix with just reads from sampleID
        sample_etpyM = [sample_etpyM; etpyM(index,:)];
    end
    % add column to start counting
    countStates = zeros(1,size(states,1)).';
    
    % go through the sample's matrix and count the observed states
    for r = 1:size(sample_etpyM,1)
        read = sample_etpyM(r,1:size(sample_etpyM,2)-1);
        for s = 1:size(states,1)
            if read == states(s,:)
                countStates(s) = countStates(s)+1;
                break;
            end
        end
    end
    
    
    figure
    bar(countStates)
    xlabel([])
    ylabel('Occurrences','FontSize',14)
    title(['Occurrences of Each State for ',num2str(age),'-Year-Old Sample'],'FontSize',14)
    box off
    xlim([0 16])
    set(gca,'XTick',[])
    % normalize the vector of counts
    s = sum(countStates);
    prob = countStates/s;
    
    % calculate entropy with probability vector
    indiv_etpy = log2(prob).*prob;
    indiv_etpy(isnan(indiv_etpy)) = [];
    etpy = [etpy; -sum(indiv_etpy)];
    
    % label
    str = [' Entropy = ',sprintf('%.2f ',etpy)];
    annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');

end
    
end
