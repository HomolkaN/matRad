function machine = matRad_downsampleMachineFile(machine,energies,fineSamplingFirst)
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

if nargin < 2
    energies = 1:length(machine.data);
end
if nargin < 3
    fineSamplingFirst = true;
end

downsampleFactor = 0.2; % Specifies percentage of grid spacing, the smaller the broader the grid
lowerPeakArea = 0.05; % Specifies percentage of peak area
upperPeakArea = 0.04; % Specifies percentage of peak area
depthCutOff = 1.2;
% depthCutOff = 0;
% cutOff = 0.005; %cutOff data at 0.001*max(Z)

for energyIx = energies%numel(machine.data)

    if fineSamplingFirst
        newDepths = min(machine.data(energyIx).depths):0.1:max(machine.data(energyIx).depths);
        if isstruct(machine.data(energyIx).Z)
            machine.data(energyIx).Z.profileORG = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).Z.profileORG,newDepths);
            machine.data(energyIx).Z.profileAPM = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).Z.profileAPM,newDepths);
        else
            machine.data(energyIx).Z = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).Z,newDepths);
        end
        machine.data(energyIx).weight = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).weight,newDepths);
        machine.data(energyIx).sigma1 = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).sigma1,newDepths);
        machine.data(energyIx).sigma2 = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).sigma2,newDepths);
        if isfield(machine.data(energyIx),'LET')
            machine.data(energyIx).LET = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).LET,newDepths);
        end
        machine.data(energyIx).depths = newDepths';
    end

    % Define updated grid
    numOfSamplePoints = downsampleFactor*numel(machine.data(energyIx).depths);
    depths = machine.data(energyIx).depths;
    newGrid = round(max(depths)/numOfSamplePoints,2);

    % CutOff data
%     cutOffPoint = find(machine.data(energyIx).Z<cutOff*max(machine.data(energyIx).Z));
%     if ~isempty(cutOffPoint)
%        cutOffDepth = depths(cutOffPoint(1));
%     end

    % Isolate peak area to use original grid
    [~,leftPos] = min(abs(machine.data(energyIx).depths-machine.data(energyIx).peakPos*(1-lowerPeakArea)));
    [~,rightPos] = min(abs(machine.data(energyIx).depths-machine.data(energyIx).peakPos*(1+upperPeakArea)));

    % Generate new depths vector
    newDepths = matRad_interp1(depths,depths,min(depths):newGrid:max(depths));

    % Only use new grid for everything that's not peak
    newDepths(newDepths>depths(leftPos) & newDepths<depths(rightPos)) = [];
    % Use original grid for the peak area
    newDepths = [newDepths; depths(leftPos:rightPos)];
    newDepths = sort(newDepths);
    newDepths = smooth(newDepths,20);
    newDepths = round(newDepths,5);

    % Plot
%     if energyIx == 1
%         plotPoints = newDepths<machine.data(energyIx).peakPos*depthCutOff;
%         edges = 0:2:max(newDepths(plotPoints));
%         figure, yyaxis left, histogram(depths,edges), hold on, histogram(newDepths(plotPoints),edges)
%         xlim([0 max(newDepths(plotPoints))])
%         xlabel('Depth (mm)')
%         ylabel('counts')
%         yyaxis right, plot(depths,machine.data(energyIx).Z,'LineWidth',2), hold on, plot(newDepths(plotPoints),matRad_interp1(depths,machine.data(energyIx).Z,newDepths(plotPoints)),'LineWidth',2,'LineStyle','--');
%         ylabel('Dose Z (a.u.)')
%     end


    % Loop through all fields and apply downsampling
    depthSize = numel(depths);
    fnames = fieldnames(machine.data(energyIx));
    for fieldIx = 1:length(fnames)
        if isstruct(machine.data(energyIx).(fnames{fieldIx}))
            subFnames = fieldnames(machine.data(energyIx).(fnames{fieldIx}));
            for subFieldIx = 1:length(subFnames)
                if ~isstruct(machine.data(energyIx).(fnames{fieldIx}).(subFnames{subFieldIx})) && ~all(size(machine.data(energyIx).(fnames{fieldIx}).(subFnames{subFieldIx}))<depthSize)
                        machine.data(energyIx).(fnames{fieldIx}).(subFnames{subFieldIx}) = matRad_interp1(depths,machine.data(energyIx).(fnames{fieldIx}).(subFnames{subFieldIx}),newDepths,'extrap');

                        % Cutoff with depth
                        if depthCutOff>0
                            machine.data(energyIx).(fnames{fieldIx}).(subFnames{subFieldIx})(newDepths>machine.data(energyIx).peakPos*depthCutOff) = [];
                        end
                end
            end

        elseif ~all(size(machine.data(energyIx).(fnames{fieldIx}))<depthSize)
            machine.data(energyIx).(fnames{fieldIx}) = matRad_interp1(depths,machine.data(energyIx).(fnames{fieldIx}),newDepths,'extrap');
            
            % Cutoff with depth
            if depthCutOff>0
                machine.data(energyIx).(fnames{fieldIx})(newDepths>machine.data(energyIx).peakPos*depthCutOff) = [];
            end

        end
    end
    
    clear depths newDepths

end

end