% given a matrix, calculate the average methylation of each column
function avgs = avgMeth(matrix)

avgs = [];

for c = 1:size(matrix,2)
    num_0 = 0;
    num_1 = 0;
    num_NaN = 0;
    for r = 1:size(matrix,1)
        if matrix(r,c) == 0
            num_0 = num_0+1;
        elseif matrix(r,c) == 1
            num_1 = num_1+1;
        else
            num_NaN = num_NaN+1;
        end
    end
    average = num_1/(num_1+num_0);
    avgs = [avgs average];
end