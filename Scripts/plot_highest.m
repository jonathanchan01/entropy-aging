% plot average methylation for the locus with high/low values for aa/age 
% en/age correlation

correlations

%% find locus with highest aa/age en/age correlation
index = 0; % this will eventually carry the index of the locus we want

list = {'High' 'Low'};
[indx1,tf] = listdlg('PromptString','Select Avg vs. Age Correlation', ...
                     'ListString',list,'SelectionMode','single');
[indx2,tf] = listdlg('PromptString','Select Ent vs. Age Correlation', ...
                     'ListString',list,'SelectionMode','single');
if indx1==1 && indx2==1
    for i=1:size(corr,2)
        [en_value,en_index] = maxk(age_en,i);
        [aa_value,aa_index] = maxk(age_aa,i);
        if index ~= 0
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
elseif indx1==1 && indx2==2
    for i=1:size(corr,2)
        [en_value,en_index] = mink(age_en,i);
        [aa_value,aa_index] = maxk(age_aa,i);
        if index ~= 0
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
elseif indx1==2 && indx2==1
    for i=1:size(corr,2)
        [en_value,en_index] = maxk(age_en,i);
        [aa_value,aa_index] = mink(age_aa,i);
        if index ~= 0
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
elseif indx1==2 && indx2==2
    for i=1:size(corr,2)
        [en_value,en_index] = mink(age_en,i);
        [aa_value,aa_index] = mink(age_aa,i);
        if index ~= 0
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
else
    fprintf('Please choose a valid combination.')
    quit
end



fprintf('Focus on %s with avg/age correlation %d and ent/age correlation %d\n',...
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

%% plot AA and entropy with age
ages = zeros(size(metrics,2),1);
aa_values = ages;
en_values = ages;
for i=1:size(metrics,2)
    ages(i) = metrics(i).age;
    aa_values(i) = metrics(i).aa;
    en_values(i) = metrics(i).en;
end

% filter NaNs
ages2 = [];
en_values2 = [];
for i=1:size(en_values,1)
   if ~isnan(en_values(i))
       en_values2 = [en_values2; en_values(i)];
       ages2 = [ages2; ages(i)];
   end
end

figure
scatter(ages,aa_values,15,"filled")
ylabel('Average','FontSize',14)
xlabel('Age','FontSize',14)
title('Average vs. Age','FontSize',14)
p = polyfit(ages,aa_values,1);
r = p(1) * [0:90] + p(2);
R = corrcoef(ages,aa_values);
hold on
plot([0:90],r,'-');
str = [' r = ',sprintf('%.2d',R(1,2))];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
box on
hold off

figure
scatter(ages2,en_values2,15,"filled")
ylabel('Entropy','FontSize',14)
xlabel('Age','FontSize',14)
title('Entropy vs. Age','FontSize',14)
p = polyfit(ages2,en_values2,1);
r = p(1) * [0:90] + p(2);
R = corrcoef(ages2,en_values2);
hold on
plot([0:90],r,'-');
str = [' r = ',sprintf('%.2d',R(1,2))];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
box on
hold off

%% create lollipop plots for old and young sample
% isolate the data to the old and young sample
BAMtoMatrix_customSample(template,locus,'P41019061105701.bam'); % 84.4 years
BAMtoMatrix_customSample(template,locus,'P41019081706783.bam'); % 7.2 years
