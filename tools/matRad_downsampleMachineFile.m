function machine = matRad_downsampleMachineFile(machine)
% Downsampling for machine files with very long vectors that lead to memory errors
%
% call
%   machine = matRad_downsampleMachineFile(machine)
%
% input
%   machine:        matRad machine file
%
% output
%   machine:        updated matRad machine file with shorter vectors
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

downsampleFactor = 0.2; % Specifies percentage of grid spacing, the smaller the broader the grid
peakArea = 0.02; % Specifies percentage of peak area

for energyIx = 1:numel(machine.data)

    % Define updated grid
    numOfSamplePoints = downsampleFactor*numel(machine.data(energyIx).depths);
    depths = machine.data(energyIx).depths;
    newGrid = round(max(depths)/numOfSamplePoints,2);

    % Isolate peak area to use original grid
    [~,leftPos] = min(abs(machine.data(energyIx).depths-machine.data(energyIx).peakPos*(1-peakArea)));
    [~,rightPos] = min(abs(machine.data(energyIx).depths-machine.data(energyIx).peakPos*(1+peakArea)));

    % Generate new depths vector
    newDepths = matRad_interp1(depths,depths,min(depths):newGrid:max(depths));

    % Use original grid vor the peak area
    newDepths(newDepths>depths(leftPos) & newDepths<depths(rightPos)) = [];
    newDepths = [newDepths; depths(leftPos:rightPos)];
    newDepths = sort(newDepths);

    % Loop through all fields and apply downsampling
    fnames = fieldnames(machine.data(energyIx));
    for fieldIx = 1:length(fnames)
        if size(machine.data(energyIx).(fnames{fieldIx}),1) > 1

            machine.data(energyIx).(fnames{fieldIx}) = matRad_interp1(depths,machine.data(energyIx).(fnames{fieldIx}),newDepths);
        
        end
    end

    clear depths newDepths

end

end