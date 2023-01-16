function [ct, cst] = matRad_addMovement_LungPhantom(ct, cst, motionPeriod, numOfCtScen, amp,visBool)
% adds artificial sinosodal patient motion by creating a deformation vector
% field and applying it to the ct.cube by geometric transformation
%
% call
%   ct = matRad_addMovement(ct, ct.motionPeriod, ct.numOfCtScen, amp)
%
% input
%   ct:             matRad ct struct
%   cst:            matRad cst struct
%   motionPeriod:   the length of a whole breathing cycle (in seconds)
%   numOfCtScen:    number of ct phases
%   amp:            amplitude of the sinosoidal movement (in pixels)
%   visBool         boolean flag for visualization
%
%   note:           1st dim --> x LPS coordinate system
%                   2nd dim --> y LPS coordinate system
%                   3rd dim --> z LPS coordinate system
%                   a positive amplitude moves the phantom to the right,
%                   anterior, inferior
%
% output
%   ct:             modified matRad ct struct including dvf and cubes for
%                   all phases
%   cst:            modified matRad cst struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('visBool','var')
    visBool = false;
end

matRad_cfg = MatRad_Config.instance();

% book keeping
ct.motionPeriod = motionPeriod;
ct.numOfCtScen = numOfCtScen;

% set type
ct.dvfType = 'pull'; % push or pull

cubeHU_temp = cell(1,numOfCtScen);
if isfield(ct,'cube')
    cube_temp = cell(1,numOfCtScen);
end

for i = 1:numOfCtScen
    % Create new ct cube
    ct.cubeHU{i} = ct.cubeHU{1};
    if isfield(ct,'cube')
        ct.cube{i} = ct.cube{1};
    end

    % Update cst array
    for j = 1:size(cst,1)
        cst{j,4}{1,i} = cst{j,4}{1,1};
    end
end

% generate scenarios
for i = 1:numOfCtScen

    %     if isfield(ct,'hlut')
    %         padValue = min(ct.hlut(:,2));
    %     else
    %         padValue = -1024;
    %     end

    ct.dvf{i} = zeros([ct.cubeDim, 3]);

    dVec = arrayfun(@(A)  A*sin((i-1)*pi / numOfCtScen)^2, amp);

    ct.dvf{i}(:,:,:,1) = dVec(1); % deformation along x direction (i.e. 2nd coordinate in dose/ct)
    ct.dvf{i}(:,:,:,2) = dVec(2);
    ct.dvf{i}(:,:,:,3) = dVec(3);

    matRad_cfg.dispInfo('Deforming ct phase %d with [dx,dy,dz] = [%f,%f,%f] voxels\n',i,dVec(1),dVec(2),dVec(3));

    %% warp ct

    % Pad with water
    cubeHU_temp{i} = imwarp(ct.cubeHU{1}, ct.dvf{i},'FillValues',0);

    if isfield(ct,'cube')
        cube_temp{i}   = imwarp(ct.cube{1},   ct.dvf{i},'FillValues',1);
    end

    % only voxel after the lung are deformed in order to preserve the lung
    % "lungEnds" is the cutoff point of the deformation
    lungEnds = 32;
    for zi = lungEnds:ct.cubeDim(1)
        for xi = 1:ct.cubeDim(2)
            for yi = 1:ct.cubeDim(3)
                ct.cubeHU{i}(zi,xi,yi) = cubeHU_temp{i}(zi,xi,yi);
                if isfield(ct,'cube')
                    ct.cube{i}(zi,xi,yi) = cube_temp{i}(zi,xi,yi);
                end
            end
        end
    end
    

    %% warp cst
    targetIx = find(contains(cst(:,2),'target'));

    tmp = zeros(ct.cubeDim);
    tmp(cst{targetIx,4}{1}) = 1;
    tmpWarp     = imwarp(tmp, ct.dvf{i});

    cst{targetIx,4}{i} = find(tmpWarp > .5);

    %% convert dvfs to [mm]
    tmp = ct.dvf{i}(:,:,:,1);
    ct.dvf{i}(:,:,:,1) = -ct.dvf{i}(:,:,:,2) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = -tmp * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = -ct.dvf{i}(:,:,:,3) * ct.resolution.z;

    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
end



if visBool
    slice = round(ct.cubeDim(3)/2);
    figure,
    for i = 1:numOfCtScen
        clf,
        imagesc(ct.cubeHU{i}(:,:,slice))
        pause(.5);
    end
end

end




