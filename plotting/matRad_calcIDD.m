function IDD = matRad_calcIDD(doseCube,calcProfile,direction)
% function to calculate integrated depth dose (IDD) of a simple boxphantom
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

doseCube(isnan(doseCube))=0;

switch direction
    case 'y'
        if ~calcProfile
            IDD = sum(sum(doseCube,2),3);
        else
            IDD = doseCube(:,round(size(doseCube,2)/2),round(size(doseCube,3)/2));
        end
    case 'x'
        if ~calcProfile
            IDD = sum(sum(doseCube,1),3);
        else
            IDD = doseCube(round(size(doseCube,1)/2),:,round(size(doseCube,3)/2));    
        end
        IDD = permute(IDD,[2,3,1]);
    case 'z'
        if ~calcProfile
            IDD = sum(sum(doseCube,2),1);
        else
            IDD = doseCube(round(size(doseCube,1)/2),round(size(doseCube,2)/2),:);    
        end
        IDD = permute(IDD,[3,1,2]);
    otherwise
        matRad_cfg.dispError('Please choose valid direction');
end

end


