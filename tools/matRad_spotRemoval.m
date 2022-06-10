function [dij,stf] = matRad_spotRemoval(dij,w,varargin)
% matRad spot removal tool
%
% call
%   dij =           matRad_spotRemoval(dij,w)
%   [dij,stf] =     matRad_spotRemoval(dij,w,stf)
%
% Example full call for protons
%   [dij2,stf2] = matRad_spotRemoval(dij,weights,stf)
%
% input
%   dij:                    old matRad dij struct
%   w:                      optimized matRad weights
%   varargin (optional):    stf: matRad steering struct
%                           thres: threshold for weights
%
% output
%   dij:                new matRad dij struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% handle inputs
if ~isempty(varargin)
    for i = 1:nargin-2
        if isstruct(varargin{i})
            stf = varargin{i};
        elseif isscalar(varargin{i})
            thres = varargin{i};
        end
    end
end
% toggle rewriting stf if stf is an input
if exist('stf','var') && nargout > 1
    calcStf = true;
end

%% calculate new spots
% For fixed threshold mode, set threshold for spot removal to 3% of the mean weight. Testing has shown differences to
% become apparent at about 5%, therefore 3% was chosen.
if ~exist('thres','var')
    thres = 0.03;
%     thres = 0.00001;
end

% save spots that have larger weight than the set threshold
newSpots = w>thres*mean(w);
% newSpots = w>thres*max(w);


%% rewrite dij and stf with new spots
if ((sum(newSpots) ~= numel(w)) && sum(newSpots) ~= dij.totalNumOfBixels) && any(size(w)>1)
    % save new weights
    dij.cutWeights = w(newSpots);
    
    % update bixel book-keeping 
    dij.bixelNum = dij.bixelNum(newSpots);
    dij.rayNum = dij.rayNum(newSpots);
    dij.beamNum = dij.beamNum(newSpots);
    dij.totalNumOfBixels = sum(newSpots);
    
    % cut out columns in already calculated sparse matrices
    dij.physicalDose{1} = dij.physicalDose{1}(:,newSpots);
    if isfield(dij,'mAlphaDose')
        dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,newSpots);
        dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,newSpots);
    end
    if isfield(dij,'mLETDose')
        dij.mLETDose{1} = dij.mLETDose{1}(:,newSpots);
    end

    % save starting and ending indices of each beam in the weight vector
    [~,beamNumIdx] = unique(dij.beamNum);
    beamNumIdx = [beamNumIdx-1;dij.totalNumOfBixels];
    
    % loop through 
    for b = 1:dij.numOfBeams
        % calculate rays and indices in current beam
        currRaysInBeam = dij.rayNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        currBixelsInRay = dij.bixelNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        [rayCnt,rayIdx] = unique(currRaysInBeam);

        % save number of rays in current beam
        dij.numOfRaysPerBeam(b) = numel(rayCnt);

        % write new stf
        if calcStf
            % calculate number of bixels in each ray
            numOfBixelsPerRay = groupcounts(currRaysInBeam);
            stf(b).numOfBixelsPerRay = numOfBixelsPerRay';
            stf(b).totalNumOfBixels = sum(stf(b).numOfBixelsPerRay);     

            % check if any rays have completely been removed and rewrite to stf
            cutRays = ismember((1:dij.numOfRaysPerBeam(b))',rayCnt);

            if any(~cutRays)
                stf(b).ray = stf(b).ray(cutRays);
                stf(b).numOfRays = sum(cutRays);
            end

            % loop through rays and write beam parameters
            for i = 1:stf(b).numOfRays
                bixelCurrRay = currBixelsInRay(rayIdx(i):rayIdx(i)+numOfBixelsPerRay(i)-1);

                stf(b).ray(i).energy = stf(b).ray(i).energy(bixelCurrRay);
                stf(b).ray(i).focusIx = stf(b).ray(i).focusIx(bixelCurrRay);
                stf(b).ray(i).rangeShifter = stf(b).ray(i).rangeShifter(bixelCurrRay);
            end
        end
    end
    
    % update total number of rays
    dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);

    % save number of removed spots and output to console (as warning to be visible)
    dij.numOfRemovedSpots = sum(~newSpots);
    matRad_cfg.dispWarning([num2str(sum(~newSpots)),'/',num2str(numel(newSpots)) ,' spots have been removed below ',num2str(100*thres),'% of the mean weight.\n'])
else
    % output warning to console
    matRad_cfg.dispWarning('no spots have been removed.')
    dij.cutWeights = w;
end
end

