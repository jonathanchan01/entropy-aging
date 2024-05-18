% Find the correlations of each metric with age.

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/results';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/results';
fileInfo = dir(fullfile(myDir, 'chr*'));
corr = struct('locus',{},'age_ch',{},'age_aa',{},'age_en',{},'ch_aa',{},'ch_en',{},'en_aa',{});

age_ch = NaN(1,size(fileInfo,1))';
age_aa = NaN(1,size(fileInfo,1))';
age_en = NaN(1,size(fileInfo,1))';
ch_aa = NaN(1,size(fileInfo,1))';
ch_en = NaN(1,size(fileInfo,1))';
en_aa = NaN(1,size(fileInfo,1))';

for n=1:size(fileInfo,1)
   fileName = fileInfo(n).name;
   corr(n).locus = fileName;
   data = importdata(fileName);

   if size(data,2) == 4
       age = data(:,1);
       ch = data(:,2);
       aa = data(:,3);
       en = data(:,4);
       temp = corrcoef(age,ch,'rows','complete');
       corr(n).age_ch = temp(1,2);
       age_ch(n) = temp(1,2);
       temp = corrcoef(age,aa,'rows','complete');
       corr(n).age_aa = temp(1,2);
       age_aa(n) = temp(1,2);
       temp = corrcoef(age,en,'rows','complete');
       corr(n).age_en = temp(1,2);
       age_en(n) = temp(1,2);
       temp = corrcoef(ch,aa,'rows','complete');
       corr(n).ch_aa = temp(1,2);
       ch_aa(n) = temp(1,2);
       temp = corrcoef(ch,en,'rows','complete');
       corr(n).ch_en = temp(1,2);
       ch_en(n) = temp(1,2);
       temp = corrcoef(en,aa,'rows','complete');
       corr(n).en_aa = temp(1,2);
       en_aa(n) = temp(1,2);
   end
end