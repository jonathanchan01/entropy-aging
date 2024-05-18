% Focus on loci that are off-diagonal. Plot age vs. entropy and age vs.
% CHALM

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics';

correlations
%%
offdiags = struct('locus',{});
j = 1;

for i=1:size(corr,2)
    if abs(corr(i).age_en)-abs(corr(i).age_ch) > 0.4
       offdiags(j).locus = corr(i).locus;
       j = j+1;
    end
end
%%
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/results';

for i=1:size(offdiags,2)
   data_unfiltered = importdata(offdiags(i).locus);
   
   % filter all rows with any NaN
   data = [];
   if sum(sum((isnan(data_unfiltered)))) > 0
       for j=1:size(data_unfiltered,1)
           if sum(isnan(data_unfiltered(j,:))) == 0
               data = [data; data_unfiltered(j,:)];
           end
       end
   end
   
   age = data(:,1);
   ch = data(:,2);
   aa = data(:,3);
   en = data(:,4);
   
   % Age vs. CHALM scatterplot
   figure
   scatter(age,ch,15,"filled")
   xlabel('Age')
   ylabel('CHALM')
   title('Age vs. CHALM',offdiags(i).locus)
   p = polyfit(age,ch,1);
   r = p(1) * [0:90] + p(2);
   R = corrcoef(age,ch);
   R_sq = R(1,2)^2;
   R2s(n) = R_sq; 
   hold on
   plot([0:90],r,'-');
   str = [' R^2 = ',sprintf('%.2d',R_sq)];
   annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
   hold off
   
   % Age vs. Entropy scatterplot
   figure
   scatter(age,en,15,"filled")
   xlabel('Age')
   ylabel('Entropy')
   title('Age vs. Entropy',offdiags(i).locus)
   p = polyfit(age,en,1);
   r = p(1) * [0:90] + p(2);
   R = corrcoef(age,en);
   R_sq = R(1,2)^2;
   R2s(n) = R_sq; 
   hold on
   plot([0:90],r,'-');
   str = [' R^2 = ',sprintf('%.2d',R_sq)];
   annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
   hold off
   
   % Age vs. AA scatterplot
   figure
   scatter(age,aa,15,"filled")
   xlabel('Age')
   ylabel('AA')
   title('Age vs. AA',offdiags(i).locus)
   p = polyfit(age,aa,1);
   r = p(1) * [0:90] + p(2);
   R = corrcoef(age,aa);
   R_sq = R(1,2)^2;
   R2s(n) = R_sq; 
   hold on
   plot([0:90],r,'-');
   str = [' R^2 = ',sprintf('%.2d',R_sq)];
   annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');
   hold off
end
