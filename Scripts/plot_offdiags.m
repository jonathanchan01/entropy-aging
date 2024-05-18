% plot average methylation for the locus that has lowest aa/age and highest
% en/age correlation

correlations

% find locus with highest aa/age en/age correlation
index = 0; % this will eventually carry the index of the locus we want

for i=1:size(corr,2)
    [en_value,en_index] = maxk(age_en,i);
    [aa_value,aa_index] = mink(abs(age_aa),i);
    if (index ~= 0)
        break
    end
    for j=1:size(en_index,1)
        % find the first match
        if(size(find(aa_index == en_index(j)),1) ~= 0)
            index = en_index(j);
            break
        end
    end
end

fprintf('Focus on %s with age/aa correlation %d and en/aage correlation %d',...
        corr(index).locus, corr(index).age_aa, corr(index).age_en)
 
%%
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics';
load('t.mat');

loci_template = strings(size(t.loci,1),2);

for i=1:size(t.loci,1)
    loci_template(i,1) = [string(t.loci{i,:})];
    loci_template(i,2) = [string(t.sequences{i,:})];
end

index = find(loci_template(:,1) == corr(index).locus);
locus = char(loci_template(index,1));
template = char(loci_template(index,2));

% get AA, CH, and EN metrics
metrics = BAMtoMatrix(template,locus);

%% plot AA aand entropy with age
ages = zeros(size(metrics,2),1);
aa_values = ages;
en_values = ages;
for i=1:size(metrics,2)
    ages(i) = metrics(i).age;
    aa_values(i) = metrics(i).aa;
    en_values(i) = metrics(i).en;
end

% filter out NaNs
ages1 = [];
aa_values2 = [];
ages2 = [];
en_values2 = [];
for i=1:size(en_values,1)
   if ~isnan(en_values(i))
       en_values2 = [en_values2; en_values(i)];
       ages2 = [ages2; ages(i)];
   end
   
   if ~isnan(aa_values(i))
       aa_values2 = [aa_values2; aa_values(i)];
       ages1 = [ages1; ages(i)];
   end
end

figure
scatter(ages1,aa_values2,15,"filled")
ylabel('AA')
xlabel('age')
title('AA vs. age for off-diagonal locus')
p = polyfit(ages1,aa_values2,1);
r = p(1) * [0:90] + p(2);
R = corrcoef(ages1,aa_values2);
R_sq = R(1,2)^2;
R2s(n) = R_sq; 
hold on
plot([0:90],r,'-');
str = [' R^2 = ',sprintf('%.2d',R_sq)];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
hold off

figure
scatter(ages2,en_values2,15,"filled")
ylabel('Entropy')
xlabel('age')
title('Entropy vs. age for off-diagonal locus')
p = polyfit(ages2,en_values2,1);
r = p(1) * [0:90] + p(2);
R = corrcoef(ages2,en_values2);
R_sq = R(1,2)^2;
R2s(n) = R_sq; 
hold on
plot([0:90],r,'-');
str = [' R^2 = ',sprintf('%.2d',R_sq)];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
hold off

%% create lollipop plots for old and young sample
% isolate the data to the old and young sample (age 7.2 and 78.2, the
% latter with AA .3)
BAMtoMatrix_customSample(template,locus,'P41019081703465.bam');
BAMtoMatrix_customSample(template,locus,'P41019081706783.bam');
