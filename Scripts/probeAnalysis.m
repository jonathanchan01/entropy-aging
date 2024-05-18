% For all the template reads in templates.txt, extract the three metrics:
% CHALM, Average of Averages, entropy
cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data';
myDir = '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data';
addpath '/Users/jonathanchan/Documents/MATLAB/epigenetics/Scripts'

% get starting information: templates, loci, and age
load('t.mat');

% forward or reverse strand analysis?
list1 = {'forward strand','reverse strand'};
[indx1,tf] = listdlg('ListString',list1,'SelectionMode','single');
if indx1 == 1
    reverse = false;
elseif indx1 == 2
    reverse = true;
end

% use each template and loci pairing to generate metrics
results = struct('locus',{},'ch',{},'aa',{},'en',{},'age',{});
for n = 1:size(t.sequences,1)
    cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Scripts'
    template = t.sequences{n};
    locus = t.loci{n};
    
    % no methylation data
    if size(find(template == 'G')) == [1,0]
       fprintf(2,'\nNo cpg sites at locus %s \n', locus);
       for r = 1:21
           results(r).ch = NaN;
           results(r).aa = NaN;
           results(r).en = NaN;
           results(r).age = NaN;
       end
    else % there is data
        metrics = BAMtoMatrix(template,locus,reverse);

        % store results
        results(n).locus=locus;
        
        for l = 1:size(metrics,2)
           results(l).locus = locus;
        end
        
        for c = 1:size(metrics,2)
            results(c).ch = metrics(c).ch;
        end
        
        for a = 1:size(metrics,2)
            results(a).aa = metrics(a).aa;
        end
        
        for e = 1:size(metrics,2)
            results(e).en = metrics(e).en;
        end
        
        for ag = 1:size(metrics,2)
            results(ag).age = metrics(ag).age;
        end
    end
    
    fprintf(1,'That was locus number %i.\n', n);

    % write out results as code runs
    cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Data/addtl_results'
    
    A = [results.age; results.ch; results.aa; results.en];    
    
    % each column is a metric
    fid = fopen(locus, 'w');
    fprintf(fid,'%12.4f %12.4f %12.4f %12.4f\r\n',A);
    fclose(fid);
    
end

cd '/Users/jonathanchan/Documents/MATLAB/epigenetics/Scripts'

