function resultGUI = matRad_recalcRBEfromLET(resultGUI,modelName,alphaX,betaX)

matRad_cfg = MatRad_Config.instance();

if nargin <=2
    alphaX = 0.1;
    betaX = 0.05;
end

if isfield(resultGUI,'LET')

else
    numOfBeams = sum(contains(fieldnames(resultGUI),'LET_beam'));
end

% if isfield(resultGUI,'w')
%     wBeam = resultGUI.w;
% else
%     wBeam = ones(resultGUI.totalNumOfBixels);
% end
if isfield(resultGUI,'LET')
    numOfVoxel = numel(resultGUI.LET);
    numOfBeams = sum(contains(fieldnames(resultGUI),'LET_beam'));
    doseGrid = size(resultGUI.LET);
    dij = 0;
elseif isfield(resultGUI,'mLETDose')
    numOfVoxel = size(resultGUI.mLETDose{1},1);
    numOfBeams = resultGUI.numOfBeams;
    numOfBeamlets = size(resultGUI.mLETDose{1},2);
    doseGrid = resultGUI.doseGrid.dimensions;
    dij = 1;
end

if ~isscalar(alphaX) && ~isscalar(betaX)
    ax = alphaX;
    bx = betaX;

    alphaX = unique(ax);
    alphaX = alphaX(alphaX ~= 0);

    betaX = unique(bx);
    betaX = betaX(betaX ~= 0);

    if length(alphaX) > 1
        matRad_cfg.dispError('Multiple alpha/beta ratios in one cube is not yet implemented');
    end
else
    ax = alphaX * ones(prod(doseGrid),1);
    bx = betaX * ones(prod(doseGrid),1);
end

ABratio = alphaX ./ betaX;

modelParam.p0_MCN   = 0.999064;     % according to https://www.ncbi.nlm.nih.gov/pubmed/26459756
modelParam.p1_MCN   = 0.35605;
modelParam.p2_MCN   = 1.1012;
modelParam.p3_MCN   = -0.0038703;

modelParam.p0_WED   = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
modelParam.p1_WED   = 0.434;
modelParam.p2_WED   = 1;

if ~dij
    mAlphaDose = spalloc(numOfVoxel,numOfBeams,1);
    mSqrtBetaDose = spalloc(numOfVoxel,numOfBeams,1);

    for i = 1:numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = ([1:numOfBeams]' == i);
    end
    beamInfo(numOfBeams+1).suffix = '';
    beamInfo(numOfBeams+1).logIx  = true(numOfBeams,1);

    for i = 1:length(beamInfo)-1
        LET = resultGUI.(['LET', beamInfo(i).suffix]);
        LET(isnan(LET)) = 0;

        % consider biological optimization
        if exist('bx','var')
            ix = bx~=0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;
        else
            ix = ones(1,prod(doseGrid));
        end

        switch modelName
            case 'MCN'

                RBEmax     = modelParam.p0_MCN + ((modelParam.p1_MCN * LET )./ ABratio);
                RBEmin     = modelParam.p2_MCN + (modelParam.p3_MCN  * sqrt(ABratio) .* LET);

                resultGUI.(['alpha_recalc_MCN', beamInfo(i).suffix]) = zeros(doseGrid);
                resultGUI.(['beta_recalc_MCN', beamInfo(i).suffix]) = zeros(doseGrid);

                resultGUI.(['alpha_recalc_MCN', beamInfo(i).suffix])(ix) = RBEmax(ix)    .* alphaX;
                resultGUI.(['beta_recalc_MCN', beamInfo(i).suffix])(ix)  = RBEmin(ix).^2 .* betaX;

            case 'WED'

                RBEmax     = modelParam.p0_WED + ((modelParam.p1_WED * LET )./ ABratio);
                RBEmin     = modelParam.p2_WED;

                resultGUI.(['alpha_recalc_WED', beamInfo(i).suffix])(ix) = RBEmax(ix)    .* alphaX;
                resultGUI.(['beta_recalc_WED', beamInfo(i).suffix])(ix)  = RBEmin(ix).^2 .* betaX;

        end

        mAlphaDose(:,i)         = reshape(resultGUI.(['alpha_recalc_' modelName, beamInfo(i).suffix]) .* resultGUI.(['physicalDose', beamInfo(i).suffix]),[],1);
        mSqrtBetaDose(:,i)      = reshape(sqrt(resultGUI.(['beta_recalc_' modelName, beamInfo(i).suffix])) .* resultGUI.(['physicalDose', beamInfo(i).suffix]),[],1);

    end

    for i = 1:length(beamInfo)

        % consider biological optimization
        if exist('bx','var')
            ix = bx~=0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;
        else
            ix = ones(1,prod(doseGrid));
        end

        resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])       = full(mAlphaDose * beamInfo(i).logIx + (mSqrtBetaDose * beamInfo(i).logIx).^2);
        resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])       = reshape(resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix]),doseGrid);

        resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])        = zeros(size(resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])));
        resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])(ix)    = (sqrt(ax(ix).^2 + 4 .* bx(ix) .* resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])(ix)) - ax(ix))./(2.*bx(ix));

        resultGUI.(['RBE_recalc_' modelName, beamInfo(i).suffix])          = resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);
    end

else

    LET = resultGUI.mLETDose{1} ./ resultGUI.physicalDose{1};
    LET(isnan(LET)) = 0;
    LET = full(LET);

    resultGUI.(['alpha_recalc_' modelName]){1} = sparse(zeros(size(resultGUI.physicalDose{1})));
    resultGUI.(['beta_recalc_' modelName]){1} = sparse(zeros(size(resultGUI.physicalDose{1})));

    for beamlet = 1:numOfBeamlets
        
        % consider biological optimization
        if exist('bx','var')
            ix = bx~=0 & resultGUI.physicalDose{1}(:,beamlet) > 0;
            ix = full(ix);
        else
            ix = ones(1,prod(doseGrid));
        end

        switch modelName
            case 'MCN'

                RBEmax     = modelParam.p0_MCN + ((modelParam.p1_MCN * LET(:,beamlet) )./ ABratio);
                RBEmin     = modelParam.p2_MCN + (modelParam.p3_MCN  * sqrt(ABratio) .* LET(:,beamlet));

                resultGUI.(['alpha_recalc_' modelName]){1}(ix,beamlet) = RBEmax(ix)    .* alphaX;
                resultGUI.(['beta_recalc_' modelName]){1}(ix,beamlet)  = RBEmin(ix).^2 .* betaX;

            case 'WED'

                RBEmax     = modelParam.p0_WED + ((modelParam.p1_WED * LET(:,beamlet) )./ ABratio);
                RBEmin     = modelParam.p2_WED;

                resultGUI.(['alpha_recalc_' modelName]){1}(ix,beamlet) = RBEmax(ix)    .* alphaX;
                resultGUI.(['beta_recalc_' modelName]){1}(ix,beamlet)  = RBEmin(ix).^2 .* betaX;

        end
    end

    resultGUI.(['mAlphaDose_recalc_' modelName]){1}      = resultGUI.(['alpha_recalc_' modelName]){1} .* resultGUI.physicalDose{1};
    resultGUI.(['mSqrtBetaDose_recalc_' modelName]){1}   = sqrt(resultGUI.(['alpha_recalc_' modelName]){1}) .* resultGUI.physicalDose{1};

    resultGUI = rmfield(resultGUI,{['alpha_recalc_' modelName],['beta_recalc_' modelName]});
end

if isfield(resultGUI,'RBE_model')
    if iscell(resultGUI.RBE_model)
        resultGUI.RBE_model{end+1} = ['recalc_' modelName];
    else
        resultGUI.RBE_model = cellstr(['recalc_' modelName]);
    end
else
    resultGUI.RBE_model{1} = ['recalc_' modelName];
end

resultGUI = orderfields(resultGUI);
end