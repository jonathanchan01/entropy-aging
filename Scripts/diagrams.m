% Make histograms for those correlations. Do the same for 
% cross-correlations.

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics';
correlations

edges = [-1:0.1:1];

% histogram to visualize distribution of correlation coefficients
figure
fig1 = tiledlayout(3,3);
nexttile
histogram(age_ch,edges)
xlabel('Correlations','FontSize',12)
ylabel('Counts','FontSize',12)
title('Distribution of Correlations: Age vs. CHALM','FontSize',12)
box on
ylim([0 350])
xlim([-1 1])
box(gca,'off')

nexttile
histogram(age_aa,edges)
xlabel('Correlations','FontSize',12)
ylabel('Counts','FontSize',12)
title('Distribution of Correlations: Age vs. Average','FontSize',12)
box on
ylim([0 350])
xlim([-1 1])
box(gca,'off')

nexttile
histogram(age_en,edges)
xlabel('Correlations','FontSize',12)
ylabel('Counts','FontSize',12)
title('Distribution of Correlations: Age vs. Entropy','FontSize',12)
box on
ylim([0 350])
xlim([-1 1])
box(gca,'off')

%% scatterplots to compare r-values
nexttile
scatter(age_aa,age_ch,15,"filled")
refline(1,0)
xlabel('Age vs. Average','FontSize',12)
ylabel('Age vs. CHALM','FontSize',12)
title('Average and CHALM method comparison','FontSize',12)

nexttile
scatter(age_aa,age_en,15,"filled")
refline(1,0)
xlabel('Age vs. Average','FontSize',12)
ylabel('Age vs. Entropy','FontSize',12)
title('Average and Entropy method comparison','FontSize',12)

nexttile
scatter(age_ch,age_en,15,"filled")
refline(1,0)
xlabel('Age vs. CHALM','FontSize',12)
ylabel('Age vs. Entropy','FontSize',12)
title('CHALM and Entropy method comparison','FontSize',12)
