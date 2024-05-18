% Given the binary matrices for all samples at a specific locus, return a
% matrix for entropy calculations.

% Stack all matrices into one big matrix. Only keep the rows and columns
% that have ≥0.5 ocupancy, keeping track of which read belongs to which
% sample. Keep the 4 columns with the least amount of NaNs and the rows
% with at most one NaN. Impute the matrix. Return the matrix and a vector
% of counts of reads across each sample.

function [final_etpyM, countSampleIDs, m_size] = makeEtpyMatrix(m)

m_size = size(m,2);

% stack matrices
bigMatrix = m(1).matrix;
for n = 1:size(m,2)-1
    bigMatrix = [bigMatrix; m(n+1).matrix];
end

if size(bigMatrix) == [0,0] % no data
    final_etpyM = [];
    countSampleIDs = zeros(79,1);
    return
end

%% filter based on occupancy

% calculate vertical occupancy of matrix:
% (values)/(values + NaN) for each column
occ_vert = [];
for c = 1:size(bigMatrix,2)
   num_count = 0;
   nan_count = 0;
   for r = 1:size(bigMatrix,1)
       if bigMatrix(r,c) == 1 || bigMatrix(r,c) == 0
           num_count = num_count + 1;
       elseif isnan(bigMatrix(r,c))
           nan_count = nan_count + 1;
       end
   end
   occ_vert = [occ_vert num_count/(num_count+nan_count)];
end

% cols with ≥0.5 occupancy -> high occupancy matrix (HOM)

HOM = [];
HOSites = [];

% get cpgSites; i is the index of the first valid vector of cpgSites in m
i = 1;
while size(m(i).cpgSites) == [0,0]
    i = i+1;
end

cpgSites = m(i).cpgSites;

for c = 1:size(occ_vert,2)
   if occ_vert(1,c) >= 0.5
       HOSites = [HOSites cpgSites(c)];
       HOM = [HOM bigMatrix(:,c)];
   end
end

% if there are less than 4 cols to work with, stop
if size((find(occ_vert>=.5)),2) < 5
    final_etpyM = [];
    countSampleIDs = zeros(79,1);
    return
end

% add a column that keeps track of sampleID
sampleID = 1;
IDcol = [];
for c = 1:size(m,2)
    numCols = size(m(c).matrix,1);
    for nn = 1:numCols
        IDcol = [IDcol; sampleID];
    end
    sampleID = sampleID + 1;
end

HOM = [HOM, IDcol];

% calculate horizontal occupancy of matrix:
% (values)/(values + NaN) for data in each row
occ_horiz = [];
for r = 1:size(HOM,1)
   num_count = 0;
   nan_count = 0;
   for c = 1:size(HOM,2)-1
       if HOM(r,c) == 1 || HOM(r,c) == 0
           num_count = num_count + 1;
       elseif isnan(HOM(r,c))
           nan_count = nan_count + 1;
       end
   end
   occ_horiz = [occ_horiz; num_count/(num_count+nan_count)];
end 

% extract rows from HOM with ≥0.5 occupancy -> HOM_final
HOM_final = [];
for r = 1:size(occ_horiz,1)
   if occ_horiz(r,1) >= 0.5
       HOM_final = [HOM_final; HOM(r,:)];
   end
end

%% keep 4 highest-occupied columns and impute

% keep the 4 columns that have the least amount of NaNs
column_nans = [];
for n = 1:size(HOM_final,2)-1
    column_nans = [column_nans, sum(isnan(bigMatrix(: ,n)))];
end

% keep track of the corresponding column number
colnums = [1:size(column_nans,2)];

% remove columns from colnum_nans and colnums until there are only 4 left
s = size(column_nans,2);
for n = 4:s-1
    max = maxk(column_nans,1);
    column_nans_temp = [];
    colnums_temp = [];
    for pos = 1:size(column_nans,2)
        % reconstruct the vector just without that key value
        if max ~= column_nans(pos)
            column_nans_temp = [column_nans_temp, column_nans(pos)];
            colnums_temp = [colnums_temp, colnums(pos)];
        else
            max = -1;
        end
    end
    
    column_nans = column_nans_temp;
    colnums = colnums_temp;
end

% construct the final matrix to do entropy calculations on
etpyM = [];
for n = 1:size(colnums,2)
   etpyM = [etpyM, HOM_final(:,colnums(n))];
end

etpyM = [etpyM, HOM_final(:,size(HOM_final,2))];

% throw out rows with more than one missing value
clean_etpyM = [];
for n = 1:size(etpyM,1)
    if sum(isnan(etpyM(n,:))) <= 1
        clean_etpyM = [clean_etpyM; etpyM(n,:)];
    end
end

% make temp matrix for imputing
etpyM_imp = [];
for c = 1:size(clean_etpyM,2)-1
    etpyM_imp = [etpyM_imp, clean_etpyM(:,c)];
end

etpyM_imp = knnimpute(etpyM_imp);

% add back the label column
final_etpyM = [etpyM_imp, clean_etpyM(:,size(clean_etpyM,2))];

% count the number of occurrences in each sample
countSampleIDs = zeros(size(m,2),1);
sampleID = final_etpyM(:,size(final_etpyM,2));
for sample = 1:size(m,2)
    count = arrayfun(@(x) length(find(sampleID(:,x)==sample)),1:size(sampleID,2));
    countSampleIDs(sample) = count;
end

end
