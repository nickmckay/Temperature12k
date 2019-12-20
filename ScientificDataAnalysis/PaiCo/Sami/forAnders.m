%calculate nests for Anders
clear
%records to omit consecutively
to_omit=[27 24 46 18 17 29 50 16 5 8 6 13 9 30 37 19 28]';

to_omit_i=cell(length(to_omit),1);
to_include_i=to_omit_i;
to_omit_i{2,1}=to_omit(1);
for i=3:length(to_omit)
    if i<5
    to_omit_i{i,1}=[to_omit_i{i-1,1} to_omit(i-1)];
    elseif i==5
            to_omit_i{i,1}=[to_omit_i{i-1,1} to_omit(i-1) to_omit(i)];
    else
                    to_omit_i{i,1}=[to_omit_i{i-1,1} to_omit(i)];
    end
    
    
end
load OrigRecon
load ~/Dropbox/ArcPaicoData.mat

for i=1:length(to_include_i)
    to_include_i{i,1}=setdiff(1:length(arc.proxy),to_omit_i{i,1});
end

%%

arc.target.times = 00:2000;
arc.instrumental.noisestd = .2224;
arc.instrumental.data = OrigRecon.Arc.Temp(1851:end)';
arc.instrumental.times = OrigRecon.Arc.Year(1851:end)';

options.resample.number = 0;
options.resample.parallelize = true;        
options.filecache = false;  
options.damping = [];
options.errorTolerance = 1e-8;
options.maxIters = 1e2;
options.heuristicstart = true;
options.regcov = 1000;
delete(gcp)
matlabpool local 4
addpath /home/nick/Dropbox/ML_Scripts/Sami/common
pctRunOnAll javaaddpath /home/nick/Dropbox/ML_Scripts/Sami/common
tracker();

for i=1:length(to_omit_i);
arcnew=arc;
arcnew.proxy=arc.proxy(to_include_i{i,1});
result = paico(arcnew,options,@tracker);

nest(:,i)=result.signal';
i
end