classdef matRad_MCsquareBaseData < matRad_MCemittanceBaseData
    % matRad_MCsquareBaseData class for calculating MCsquare base data and
    % writing it to a .txt file, for MCsquare to use
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
    
    methods (Access = public)
        function obj = matRad_MCsquareBaseData(machine,stf)
            %Call matRad_MCemmitanceBaseData constructor
            if nargin < 2
                stf = [];
            end
            
            obj = obj@matRad_MCemittanceBaseData(machine, stf);
        end
                
        function obj = writeMCsquareData(obj,filepath)
            %function that writes a data file containing Monte Carlo base
            %data for a simulation with MCsquare
                       
            %write MCsqaure data base file
            try
                
                fileID = fopen(filepath,'w');
                
                %Header
                %fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
                fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
                fprintf(fileID,'# %s\n', obj.machine.meta.description);
                fprintf(fileID,'# created by %s on %s\n\n', obj.machine.meta.created_by, obj.machine.meta.created_on);
                
                fprintf(fileID,'Nozzle exit to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.nozzleToIso);
                
                fprintf(fileID,'SMX to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smx);
                
                fprintf(fileID,'SMY to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smy);
                
                for i = 1:length(obj.rangeShifters)
                    raShi = obj.rangeShifters(i);
                    fprintf(fileID,'Range Shifter parameters\n');
                    fprintf(fileID,'RS_ID = %d\n',raShi.ID);
                    fprintf(fileID,'RS_type = binary\n');
                    fprintf(fileID,'RS_material = 64\n'); %Material ID Hardcoded for now (PMMA)
                    fprintf(fileID,'RS_density = 1.19\n'); %Maetiral density Hardcoded for now (PMMA)
                    fprintf(fileID,'RS_WET = %f\n\n',raShi.eqThickness);
                end
                    
                % Get unique energies with focusIndices
                uniqueEnergyFocus = unique(obj.focusTable,'rows');

                % Write number of unique energies
                fprintf(fileID,'Beam parameters\n%d energies\n\n',size(uniqueEnergyFocus,1));
                
                % Write header fields
                fn = fieldnames(obj.monteCarloData)';
                fprintf(fileID, '%s\t', fn{:});
                fprintf(fileID, '\n');

                % Write energies
                for k = 1:size(uniqueEnergyFocus,1)
                    currentEnergyFocus = table2array(uniqueEnergyFocus(k,:));
                    [~,energyIx] = intersect(unique([obj.monteCarloData.NominalEnergy]),currentEnergyFocus(1));

                    fprintf(fileID, '%g\t', structfun(@(fld) fld(currentEnergyFocus(2)), obj.monteCarloData(energyIx))');
                    fprintf(fileID, '\n');
                end
                
                fclose(fileID);
                
                obj.bdl_path = filepath;
                
            catch MException
                error(MException.message);
            end
        end        
    end    
end

