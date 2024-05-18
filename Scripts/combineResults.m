% combine results. For each file in results, find the corresponding file in
% addtl_results and append the values of addtl_results

%% get addtl data
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/addtl_results';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/addtl_results';
addtl_fileInfo = dir(fullfile(myDir, 'chr*'));
n = size(addtl_fileInfo,1);

addtl_data = struct('locus',{},'age',{},'aa',{},'ch',{},'en',{});
for k=1:n
    addtl_data(k).locus = addtl_fileInfo(k).name;
    tempdata = importdata(addtl_fileInfo(k).name);
    datasize = size(tempdata,2);
    switch datasize
        case 4 % every metric is there
            addtl_data(k).age = tempdata(:,1);
            addtl_data(k).ch = tempdata(:,2);
            addtl_data(k).aa = tempdata(:,3);
            addtl_data(k).en = tempdata(:,4);
        case 3 % entropy is all NaN
            y = tempdata(:,1);
            addtl_data(k).age = tempdata(:,1);
            addtl_data(k).ch = tempdata(:,2);
            addtl_data(k).aa = tempdata(:,3);
        case 2 % entropy and CHALM NaN
            y = tempdata(:,1);
            addtl_data(k).age = tempdata(:,1);
            addtl_data(k).ch = tempdata(:,2);
        case 0
            addtl_data(k).age = nan(1,size(tempdata,1))';
            addtl_data(k).ch = nan(1,size(tempdata,1))';
            addtl_data(k).aa = nan(1,size(tempdata,1))';
            addtl_data(k).en = nan(1,size(tempdata,1))';
    end
end

%% get original data

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/results';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/results';
fileInfo = dir(fullfile(myDir, 'chr*'));
n = size(fileInfo,1);

data = struct('locus',{},'age',{},'aa',{},'ch',{},'en',{});
for k=1:n
    data(k).locus = fileInfo(k).name;
    tempdata = importdata(fileInfo(k).name);
    datasize = size(tempdata,2);
    switch datasize
        case 4 % every metric is there
            data(k).age = tempdata(:,1);
            data(k).ch = tempdata(:,2);
            data(k).aa = tempdata(:,3);
            data(k).en = tempdata(:,4);
        case 3 % entropy is all NaN
            y = tempdata(:,1);
            data(k).age = tempdata(:,1);
            data(k).ch = tempdata(:,2);
            data(k).aa = tempdata(:,3);
        case 2 % entropy and CHALM NaN
            y = tempdata(:,1);
            data(k).age = tempdata(:,1);
            data(k).ch = tempdata(:,2);
        case 0
            data(k).age = nan(1,size(tempdata,1))';
            data(k).ch = nan(1,size(tempdata,1))';
            data(k).aa = nan(1,size(tempdata,1))';
            data(k).en = nan(1,size(tempdata,1))';
    end
end

%% combine
for k=1:n
    if ~isempty(data(k).age) && ~isempty(addtl_data(k).age)
        data(k).age(80:100) = addtl_data(k).age;
    end
    if ~isempty(data(k).aa) && ~isempty(addtl_data(k).aa)
        data(k).aa(80:100) = addtl_data(k).aa;
    end
    if ~isempty(data(k).ch) && ~isempty(addtl_data(k).ch)
        data(k).ch(80:100) = addtl_data(k).ch;
    end
    if ~isempty(data(k).en) && ~isempty(addtl_data(k).en)
        data(k).en(80:100) = addtl_data(k).en;
    end
    if ~isempty(data(k).en) && ~isempty(addtl_data(k).en)
        data(k).en(80:100) = addtl_data(k).en;
    end
end