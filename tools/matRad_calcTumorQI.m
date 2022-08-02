function outTable = matRad_calcTumorQI(cst,varargin)
% volumePTV = Volume of the PTV segmentation in [mL] or [cc]
matRad_cfg = MatRad_Config.instance();
field = {'PTV','target'};
meanDose = [];

for i = 1:nargin-1
    isMat{i} = structfun(@(x) ndims(x)==3, varargin{i});
    f_orig{i} = fieldnames(varargin{i});
    f_orig{i} = f_orig{i}(isMat{i});
end
f = unique(cat(1, f_orig{:}));

idx = cst(contains(cst(:,2),field),[2,4]);
tabHeader = {'Field','mean PTV'};
tabHeader2 = [];

for i = 1:nargin-1
    clear r
    f1 = fieldnames(varargin{i});
    f1 = f1(isMat{i});
    y = ismember(f,f1);

    tabHeader2 = [tabHeader2 , {inputname(i+1)}];

    currMean = structfun(@(x)  matRad_calcMean(x,idx)  , varargin{i});
    r(y) = currMean(currMean>-100);
    r(~y) = NaN;
    meanDose = [meanDose r'];
end

outTable = table(f,meanDose);
outTable = splitvars(outTable,'meanDose','NewVariableNames',tabHeader2);
end

function meanS = matRad_calcMean(x,idx)
     if ndims(x)==3
        meanS = mean(x(idx{1,2}{1}));
     else
        meanS = -100;
     end
end