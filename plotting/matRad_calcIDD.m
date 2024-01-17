function outputDD = matRad_calcIDD(doseCube,calcProfile,direction,slice,averageProfiles)
% function to calculate depth dose (DD) of a simple boxphantom
% for beams in x-, y- or z-direction
%
% call
%   IDD = matRad_calcIDD(doseCube,direction)
%
% input
%   doseCube:       calculated dose cube
%   direction:      'x','y' or 'z', direction of the beam
%
% output
%   IDD:            depth dose vector
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isstruct(doseCube)
    error('Cube was entered as struct.')
end
if nargin < 3
    direction = 'y';
end
if nargin < 2
    calcProfile = false;
end

if nargin < 4
    slice = floor(size(doseCube)/2); %voxels
end

if nargin < 5
    averageProfiles = true;
end
doseCube(isnan(doseCube))=0;
intWidth = 2; %voxels

if calcProfile
    if averageProfiles
        % Calculate average profiles
        switch direction
            case 'y'
                outputDD = squeeze(sum(doseCube(:,slice(2)-intWidth:slice(2)+intWidth,slice(3)-intWidth:slice(3)+intWidth), [2,3]));
            case 'x'
                outputDD = squeeze(sum(doseCube(slice(1)-intWidth:slice(1)+intWidth,:,slice(3)-intWidth:slice(3)+intWidth), [1,3]));
            case 'z'
                outputDD = squeeze(sum(doseCube(slice(1)-intWidth:slice(1)+intWidth,slice(2)-intWidth:slice(2)+intWidth,:), [1,2]));
            otherwise
                matRad_cfg.dispError('Please choose valid direction');
        end

    else
        % Calculate single profiles
        switch direction
            case 'y'
                outputDD = squeeze(doseCube(:,slice(2),slice(3)));
            case 'x'
                outputDD = squeeze(doseCube(slice(1),:,slice(3)));
            case 'z'
                outputDD = squeeze(doseCube(slice(1),slice(2),:));
            otherwise
                matRad_cfg.dispError('Please choose valid direction');
        end
    end
else
    % Calculate IDD
    switch direction
        case 'y'
            outputDD = squeeze(sum(doseCube,[2,3]));
        case 'x'
            outputDD = squeeze(sum(doseCube,[1,3]));
        case 'z'
            outputDD = squeeze(sum(doseCube,[1,2]));
        otherwise
            matRad_cfg.dispError('Please choose valid direction');
    end
end


end

