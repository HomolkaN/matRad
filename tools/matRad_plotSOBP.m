function matRad_plotSOBP(ct, stf, cube1, cube2)
% matRad SOBP plotting function
%
% call
%   matRad_plotSOBP(ct, stf, cube1, cube2)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   cube1:          dose cube 1
%   cube2:          dose cube 2 of same size for comparison (optional)
%
% output
%   figure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('cube2')
    checkSize = size(cube1) == size(cube2);
    if ~all(checkSize)
        error('Dose cubes are not the same size!')
    end
end
for i = 1:length(stf)
    for k = 1:length(stf(i).ray)
        if sum(stf(i).ray(k).rayPos_bev == [0 0 0]) == 3
            middleRayIx = k;
        end
    end
    
    [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(stf(i).isoCenter,ct.resolution,stf(i).sourcePoint,stf(i).ray(middleRayIx).targetPoint,{cube1});
    
    alphaMid = (alphas(1:end-1) + alphas(2:end))/2;
    physMid = alphaMid * d12;
    values = cube1(ix);
    
    figure;
    plot(physMid,values)
    if exist('cube2')
        values2 = cube2(ix);
        hold on
        plot(physMid,values2)
        legend('Homogeneous','Heterogeneous')
    end
    title(['Beam ',num2str(i)])
    ylabel('RBE Dose [Gy]')
    xlabel('Distance to source')
    xlim([6500 6575])
end
end



