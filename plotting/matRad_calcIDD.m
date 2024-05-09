function outputDD = matRad_calcIDD(doseCube,calcProfile,direction,rayPosition,averageProfiles,intWidth)
% function to calculate depth dose (DD) of a simple boxphantom
% for beams in x-, y- or z-direction
%
% call
%   IDD = matRad_calcIDD(doseCube,calcProfile,direction)
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
    calcProfile = true;
end

if nargin < 4
    rayPosition = floor(size(doseCube)/2); %voxels
end

if nargin < 5
    averageProfiles = true;
end
doseCube(isnan(doseCube))=0;
if nargin < 6
    intWidth = 1; %voxels (equals (2n+1)² central rays)
end

if calcProfile
    if averageProfiles
        % Calculate average profiles
        switch direction
            case 'y'
                outputDD = squeeze(sum(doseCube(:,rayPosition(2)-intWidth:rayPosition(2)+intWidth,rayPosition(3)-intWidth:rayPosition(3)+intWidth), [2,3]));
            case 'x'
                outputDD = squeeze(sum(doseCube(rayPosition(1)-intWidth:rayPosition(1)+intWidth,:,rayPosition(3)-intWidth:rayPosition(3)+intWidth), [1,3]))';
            case 'z'
                outputDD = squeeze(sum(doseCube(rayPosition(1)-intWidth:rayPosition(1)+intWidth,rayPosition(2)-intWidth:rayPosition(2)+intWidth,:), [1,2]))';
            otherwise
                matRad_cfg.dispError('Please choose valid direction');
        end

    else
        % Calculate single profiles
        switch direction
            case 'y'
                outputDD = squeeze(doseCube(:,rayPosition(2),rayPosition(3)));
            case 'x'
                outputDD = squeeze(doseCube(rayPosition(1),:,rayPosition(3)));
            case 'z'
                outputDD = squeeze(doseCube(rayPosition(1),rayPosition(2),:));
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

