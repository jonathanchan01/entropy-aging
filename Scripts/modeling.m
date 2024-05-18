% Perform elastic net regression or neural network regression on a certain metric
% or combination of metrics, based on user input.

%% Use master_dataset, create the matrix X of observations and vector y of ages
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data';
load('master_dataset.mat');
n = size(data,2);

y = [data(1).age];
X = [];

% user picks which regression to perform
list1 = {'Average Methylation' 'CHALM' 'Methylation Entropy' 'All'};
[indx1,tf] = listdlg('ListString',list1,'SelectionMode','single');
list2 = {'Elastic Net','Neural Net'};
[indx2,tf] = listdlg('ListString',list2,'SelectionMode','single');


% let y be the vector of ages, X be the table of all aa/ch/en values
% rows of X are samples, columns are loci
for k=1:n
    % either entropy is missing, or everything is missing
    if ~isempty(data(k).age) % entropy is maybe missing
        if size(data(k).age,1) == 79
            continue
        end
        switch indx1
            case 1
                X = [X data(k).aa];
            case 2
                X = [X data(k).ch];
            case 3
                if size(data(k).en,1) == 79
                    continue
                else
                    X = [X data(k).en];
                end
            case 4
                if size(data(k).en,1) == 79
                    continue
                else
                    X = [X data(k).aa data(k).ch data(k).en];
                end
        end
    end
end

%% in order to use knnimpute on a column, it must have at least one entry without nan
% make a clean matrix X_clean s.t. imputation can be performed
thresh = size(y,1)-1; 
X_temp = X;

while size(find(sum(isnan(X_temp),2)==0),1) == 0 % every column has an nan
    % find columns with threshold amount of nans
    nanmatrix = isnan(X_temp);
    numnans = sum(nanmatrix);
    badcols = find(numnans > thresh);

    % X_clean is the new matrix without those columns
    X_clean = zeros(size(X_temp,1),size(X_temp,2)-size(badcols,2));    
    counter = 1;
    for j=1:size(X_temp,2)
        if counter <= size(badcols,2) && j==badcols(counter)
            counter = counter+1;
        else
            X_clean(:,j-counter+1)=X_temp(:,j);
        end
    end

    X_temp = X_clean;
    thresh = thresh-1;
    
    if thresh < 1
        break
    end
end

X_imputed = knnimpute(X_clean);

%% find hyperparameters
i = length(y);
yhat = nan(1,i)';
if indx2 == 2
    rng("default"); % for reproducability
    Mdl = fitrnet(X_imputed,y,"OptimizeHyperparameters","auto",...
           "HyperparameterOptimizationOptions",... % for reproducibility
           struct("AcquisitionFunctionName","expected-improvement-plus"));
    LayerSizes = Mdl.ModelParameters.LayerSizes;
    Lambda = Mdl.ModelParameters.Lambda;
    Activation = Mdl.ModelParameters.Activations;
    Standardize = Mdl.ModelParameters.StandardizeData;
end

%% Regression with LOOCV
model_coefficients = [];
    for j=1:i
       % Create test and training sets 
       yTrain = y;
       yTest = y(j);
       yTrain(j) = [];
       xTrain = X_imputed;
       xTest = X_imputed(j,:);
       xTrain(j,:) = [];

       if indx2 == 1 % do elastic net regression
           % Alpha = 1 --> lasso, 0 --> ridge 
           [B,FitInfo] = lasso(xTrain,yTrain,'Alpha',0.75,'CV',10);

           % the largest value of lambda s.t. MSE ≤ 1 s.e. of minimum MSE
           idxLambda1SE = FitInfo.Index1SE;
           coef = B(:,idxLambda1SE);
           model_coefficients = [model_coefficients coef];

           % corresponding intercept
           coef0 = FitInfo.Intercept(idxLambda1SE);

           % calculate predicted age with test data
           yhat(j) = xTest*coef + coef0;
       end
       
       if indx2 == 2 % do neural net regression
           temp_NN = fitrnet(xTrain,yTrain,'LayerSizes',LayerSizes,'Lambda',...
                     Lambda,'Activations',Activation,'Standardize',Standardize,...
                     'InitialStepSize','auto');
           yhat(j) = predict(temp_NN,xTest);
       end

       fprintf('\n%d completed of %d',j,i)
    end
 
%% Plot results

% scatterplot
figure
scatter(y,yhat,'filled')
box on
xlim([0,90])
ylim([0,90])
hold on
ylabel('Predicted Age','FontSize',14)
xlabel('Chronological Age','FontSize',14)
switch indx2
    case 1
        modeltype = 'Elastic Net';
    case 2
        modeltype = 'Neural Net';
end
switch indx1
    case 1
        title([modeltype,' Regression with Average Methylation'],'Fontsize',14)
    case 2
        title([modeltype,' Regression with CHALM'],'Fontsize',14)
    case 3
        title([modeltype,' Regression with Methylation Entropy'],'Fontsize',14)
    case 4
        title([modeltype,' Regression with All Three Metrics'],'Fontsize',14)
end


r = corrcoef(y,yhat,'rows','complete');
cc=sprintf('r = %1.3f',r(1,2));
T = text(max(get(gca, 'xlim')-2), max(get(gca, 'ylim')-2), cc);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'right');

% omit nans
plot_data = [y yhat];
plot_data(any(isnan(plot_data),2),:)=[];

coefs = polyfit(plot_data(:,1),plot_data(:,2),1);
plot(polyval(coefs,0:90))
bestfit=sprintf('y = %1.3fx + %1.3f',coefs(1),coefs(2));
T = text(min(get(gca, 'xlim')+2), max(get(gca, 'ylim')-2), bestfit);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

mae = 0;
for k=1:size(plot_data,1)
    mae = mae + abs(plot_data(k,1)-plot_data(k,2));
end

mae = mae/size(plot_data,1);
error=sprintf('MAE: %1.3f',mae);
T = text(min(get(gca, 'xlim')+2), max(get(gca, 'ylim')-7), error);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

hold off

figure
box off
ylim([-25 25])
xlim([0 90])
scatter(y,yhat-y)
ylabel('Error','FontSize',14)
xlabel('Chronological Age','FontSize',14)
switch indx1
    case 1
        title([modeltype,' Regression with Average Methylation'],'Fontsize',14)
    case 2
        title([modeltype,' Regression with CHALM'],'Fontsize',14)
    case 3
        title([modeltype,' Regression with Methylation Entropy'],'Fontsize',14)
    case 4
        title([modeltype,' Regression with All Three Metrics'],'Fontsize',14)
end