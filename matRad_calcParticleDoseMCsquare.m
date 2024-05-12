function dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad MCsquare Monte Carlo proton dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   calcDoseDirect:             binary switch to enable forward dose
%                               calculation
% output
%   dij:                        matRad dij struct
%
% References
%
%   https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4943377
%   http://www.openmcsquare.org/
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% Handle inputs
if nargin < 5
    calcDoseDirect = false;
end

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix with MCsquare...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% check if valid machine
if ~strcmp(pln.radiationMode,'protons')
    matRad_cfg.dispError('Wrong radiation modality. MCsquare only supports protons!');
end

% Load class variables in pln
% for calcDoseDirect, this is already done in superior function
if ~isfield(pln,'propMC') || ~isa(pln.propMC,'matRad_MCsquareConfig')
    pln = matRad_cfg.getDefaultClass(pln,'propMC','matRad_MCsquareConfig');
end

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,'propDoseCalc');

% override default parameters from external parameters if available
if isfield(pln,'propHeterogeneity') && ~isempty(pln.propHeterogeneity) && isprop(pln.propHeterogeneity,'sampling') && isfield(pln.propHeterogeneity.sampling,'numHistories') && isfield(ct,'modulated')
    pln.propMC.numHistories = pln.propHeterogeneity.sampling.numHistories;
    matRad_cfg.dispWarning('Histories have been overwritten by heterogeneity correction! New histories: %.2e',pln.propMC.numHistories)
end

env = matRad_getEnvironment();

%% check if binaries are available
%Executables for simulation
if ~pln.propMC.externalCalculation
    if ispc
        if exist('MCSquare_windows.exe','file') ~= 2
            matRad_cfg.dispError('Could not find MCsquare binary.\n');
        else
            mcSquareBinary = 'MCSquare_windows.exe';
        end
    elseif ismac
        if exist('MCsquare_mac','file') ~= 2
            matRad_cfg.dispError('Could not find MCsquare binary.\n');
        else
            mcSquareBinary = './MCsquare_mac';
        end
        %error('MCsquare binaries not available for mac OS.\n');
    elseif isunix
        if exist('MCsquare_linux','file') ~= 2
            matRad_cfg.dispError('Could not find MCsquare binary.\n');
        else
            mcSquareBinary = 'chmod a+x MCsquare_linux && ./MCsquare_linux';
        end
    end
end

%Mex interface for import of sparse matrix
if ~calcDoseDirect && ~matRad_checkMexFileExists('matRad_sparseBeamletsReaderMCsquare')
    matRad_cfg.dispWarning('Compiled sparse reader interface not found. Trying to compile it on the fly!');
    try
        matRad_compileMCsquareSparseReader();
    catch MException
        matRad_cfg.dispError('Could not find/generate mex interface for reading the sparse matrix. \nCause of error:\n%s\n Please compile it yourself.',MException.message);
    end
end

% set and change to MCsquare binary folder
currFolder = pwd;
fullfilename = mfilename('fullpath');
MCsquareFolder = [fullfilename(1:find(fullfilename==filesep,1,'last')) 'MCsquare' filesep 'bin'];

% cd to MCsquare folder (necessary for binary)
cd(MCsquareFolder);

%Check Materials
if ~exist([MCsquareFolder filesep 'Materials'],'dir') || ~exist(fullfile(MCsquareFolder,'Materials','list.dat'),'file')
    matRad_cfg.dispInfo('First call of MCsquare: unzipping Materials...');
    unzip('Materials.zip');
    matRad_cfg.dispInfo('Done!\n');
end

% Since MCsquare 1.1 only allows similar resolution in x&y, we do some
% extra checks on that before calling calcDoseInit. First, we make sure a
% dose grid resolution is set in the pln struct
if ~isfield(pln.propDoseCalc,'doseGrid')
    pln.propDoseCalc.doseGrid.resolution = ct.resolution;
end
if pln.propDoseCalc.doseGrid.resolution.x ~= pln.propDoseCalc.doseGrid.resolution.y
    pln.propDoseCalc.doseGrid.resolution.x = mean([pln.propDoseCalc.doseGrid.resolution.x pln.propDoseCalc.doseGrid.resolution.y]);
    pln.propDoseCalc.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.x;
    matRad_cfg.dispWarning('Anisotropic resolution in axial plane for dose calculation with MCsquare not possible\nUsing average x = y = %g mm\n',pln.propDoseCalc.doseGrid.resolution.x);
end

% Write patientID to foldername
if isfield(ct,'patientID')
    pln.propMC.MCrun_Directory = [pln.propMC.MCrun_Directory ct.patientID '_'];
end

% Set nested folder structure also for MCsquare
pln.propMC.MCrun_Directory = [pln.propMC.MCrun_Directory pln.radiationMode,'_',pln.machine,'_',date];

% Write sampleIdx to foldername
if isfield(ct,'sampleIdx')
    pln.propMC.MCrun_Directory = [pln.propMC.MCrun_Directory '_sample' num2str(ct.sampleIdx,'%02.f')];
end

% Check if Directory name is too long (yes this is a thing apparently)
if numel(pln.propMC.MCrun_Directory) > 75
    matRad_cfg.dispWarning('MCsquare MCrunDirectory too long, trying to short name.');
    pln.propMC.MCrun_Directory = erase(pln.propMC.MCrun_Directory,'_');
    if numel(pln.propMC.MCrun_Directory) > 75
        matRad_cfg.dispError('MCsquare MCrunDirectory too long.');
    end
end

% Set appropriate output Directory
pln.propMC.MCrun_Directory  = [pln.propMC.MCrun_Directory '/'];
pln.propMC.Output_Directory = [pln.propMC.MCrun_Directory 'MCsquareOutput/'];

% Clear previous data
try
    rmdir(pln.propMC.MCrun_Directory,'s')
catch
    matRad_cfg.dispInfo('Processed files could not be deleted, consider deleting them manually if an error occurs.')
end

if ~exist(['./' pln.propMC.MCrun_Directory],'dir')
    mkdir(pln.propMC.MCrun_Directory)
end

%%Now we can run calcDoseInit as usual
matRad_calcDoseInit;

%Issue a warning when we have more than 1 scenario
if dij.numOfScenarios ~= 1
    matRad_cfg.dispWarning('MCsquare is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
end

% prefill ordering of MCsquare bixels
dij.MCsquareCalcOrder = NaN*ones(dij.totalNumOfBixels,1);

% Explicitly setting the number of threads for MCsquare, 0 is all available
nbThreads = 0;

% set relative dose cutoff for storage in dose influence matrix, we use the
% default value for the lateral cutoff here
relDoseCutOff = 1 - matRad_cfg.propDoseCalc.defaultLateralCutOff;

% book keeping - this is necessary since pln is not used in optimization or
% matRad_calcCubes
if any(strcmp(pln.bioParam.model,{'constRBE','MCN','WED'}))
    dij.RBE = pln.bioParam.RBE;
    dij.RBE_model = pln.bioParam.model;
end

scenCount = 0;

% for MCsquare we explicitly downsample the ct to the dose grid (might not
% be necessary in future MCsquare versions with separated grids)
for s = 1:ct.numOfCtScen
    HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
        dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
end

% switch for using existing BDL file (e.g. to fit matRad basedata),
% or generate BDL file from matRad base data using MCsquareBDL
if isfield(pln,'loadExistingBDL') && ~isempty(pln.loadExistingBDL)
    % use existing BDL file
    bdFile = pln.loadExistingBDL;

else
    % fit and create BDL file using selected machine file
    bdFile = [machine.meta.machine '.txt'];

    % override base data in case of APM, it is not needed here
    machine.data = matRad_HeterogeneityConfig.overrideBaseData(machine.data);

    % Calculate MCsquare base data
    % Argument stf is optional, if given, calculation only for energies given in stf
    MCsquareBDL = matRad_MCsquareBaseData(machine,stf);

    %matRad_createMCsquareBaseDataFile(bdFile,machine,1);
    MCsquareBDL = MCsquareBDL.writeMCsquareData([MCsquareFolder filesep 'BDL' filesep bdFile]);
    try
        MCsquareBDL = MCsquareBDL.saveMatradMachine('savedMatRadMachine');
    catch ME
        % TODO
    end

end

for scenarioIx = 1:pln.multScen.totNumScen

    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(scenarioIx,:);
    end

    ctScen = pln.multScen.linearMask(scenarioIx,1);
    shiftScen = pln.multScen.linearMask(scenarioIx,2);
    rangeShiftScen = pln.multScen.linearMask(scenarioIx,3);

    if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)

        %For direct dose calculation
        totalWeights = [];

        %Count the scenarios
        scenCount = scenCount + 1;

        % We need to adjust the offset used in matRad_calcDoseInit
        mcSquareAddIsoCenterOffset = [dij.doseGrid.resolution.x/2 dij.doseGrid.resolution.y/2 dij.doseGrid.resolution.z/2] ...
            - [dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z];
        mcSquareAddIsoCenterOffset = mcSquareAddIsoCenterOffset - offset;

        % MCsquare settings
        MCsquareConfigFile = [pln.propMC.MCrun_Directory 'MCsquareConfig.txt'];

        pln.propMC.BDL_Machine_Parameter_File = ['BDL/' bdFile];

        pln.propMC.BDL_Plan_File = [pln.propMC.MCrun_Directory 'currBixels.txt'];
        pln.propMC.CT_File       = [pln.propMC.MCrun_Directory 'MC2patientCT.mhd'];

        pln.propMC.Num_Threads   = nbThreads;
        pln.propMC.RNG_Seed      = 1234;

        % turn simulation of individual beamlets
        pln.propMC.Beamlet_Mode = ~calcDoseDirect;
        % turn of writing of full dose cube
        pln.propMC.Dose_MHD_Output = calcDoseDirect;
        % turn on sparse output
        pln.propMC.Dose_Sparse_Output = ~calcDoseDirect;
        % set threshold of sparse matrix generation
        pln.propMC.Dose_Sparse_Threshold = relDoseCutOff;

        %Matrices for LET
        if pln.propDoseCalc.calcLET
            pln.propMC.LET_MHD_Output		 = calcDoseDirect;
            pln.propMC.LET_Sparse_Output	 = ~calcDoseDirect;
        end

        dij.beamNum = [];
        dij.rayNum = [];
        dij.bixelNum = [];

        for i = 1:length(stf)
            %Let's check if we have a unique or no range shifter, because MCsquare
            %only allows one range shifter type per field which can be IN or OUT
            %per spot
            raShiField = [];
            for j = 1:stf(i).numOfRays
                if isfield(stf(i).ray(j),'rangeShifter')
                    raShiField = [raShiField stf(i).ray(j).rangeShifter(:).ID];
                else
                    raShiField = [raShiField zeros(size(stf(i).ray(j).energies))];
                end
            end

            raShiField = unique(raShiField); %unique range shifter
            raShiField(raShiField == 0) = []; %no range shifter
            if numel(raShiField) > 1
                matRad_cfg.dispError('MCsquare does not support different range shifter IDs per field! Aborting.\n');
            end

            if ~isempty(raShiField)
                stfMCsquare(i).rangeShifterID = raShiField;
                stfMCsquare(i).rangeShifterType = 'binary';
            else
                stfMCsquare(i).rangeShifterID = 0;
                stfMCsquare(i).rangeShifterType = 'binary';
            end

            stfMCsquare(i).gantryAngle = mod(180-stf(i).gantryAngle,360); %Different MCsquare geometry
            stfMCsquare(i).couchAngle  = stf(i).couchAngle;
            stfMCsquare(i).isoCenter   = stf(i).isoCenter + mcSquareAddIsoCenterOffset;
            stfMCsquare(i).energies    = unique([stf(i).ray.energy]);
            stfMCsquare(i).SAD         = stf(i).SAD;

            % allocate empty target point container
            for j = 1:size(MCsquareBDL.focusTable,1)
                stfMCsquare(i).energyLayer(j).targetPoints   = [];
                stfMCsquare(i).energyLayer(j).numOfPrimaries = [];
                stfMCsquare(i).energyLayer(j).MU             = [];
                stfMCsquare(i).energyLayer(j).rayNum         = [];
                stfMCsquare(i).energyLayer(j).bixelNum       = [];
            end

            currBeam.bixelIx = cell2mat(cellfun(@(x) 1:x, num2cell(stf(i).numOfBixelsPerRay), 'UniformOutput', false));
            currBeam.rayIndices = repelem(1:numel(stf(i).numOfBixelsPerRay),stf(i).numOfBixelsPerRay);
            currBeam.weights = vertcat(stf(i).ray.weight);

            % Collect totalWeights
            if calcDoseDirect
                totalWeights(i) = sum(currBeam.weights);
            end

            % Get dij counter for current bea,
            dij.beamNum = [dij.beamNum, i * ones(1,stf(i).totalNumOfBixels)];
            dij.rayNum = [dij.rayNum, currBeam.rayIndices];
            dij.bixelNum = [dij.bixelNum, currBeam.bixelIx];

            % Collect energy layers for different focus indices (needed to write MCsquare data)
            for layerIx = 1:size(MCsquareBDL.focusTable,1)
                % Collect which bixels are in the current energy layer
                currEnergyLayer.isBixel = (MCsquareBDL.bixelIndices{i} == layerIx);

                % Collect bixels and rays
                stfMCsquare(i).energyLayer(layerIx).bixelNum = currBeam.bixelIx(currEnergyLayer.isBixel);
                stfMCsquare(i).energyLayer(layerIx).rayNum   = currBeam.rayIndices(currEnergyLayer.isBixel);

                if isempty(stfMCsquare(i).energyLayer(layerIx).bixelNum)
                    continue
                end
                % Collect Target points
                currEnergyLayer.targetPoints = cell2mat({stf(i).ray(stfMCsquare(i).energyLayer(layerIx).rayNum').rayPos_bev}');
                stfMCsquare(i).energyLayer(layerIx).targetPoints = [-currEnergyLayer.targetPoints(:,1) currEnergyLayer.targetPoints(:,3)];

                % Get number of primaries and MU (equal)
                stfMCsquare(i).energyLayer(layerIx).numOfPrimaries = pln.propMC.numHistories * currBeam.weights(currEnergyLayer.isBixel)';
                stfMCsquare(i).energyLayer(layerIx).MU = pln.propMC.numHistories * currBeam.weights(currEnergyLayer.isBixel)';

                %Now add the range shifter
                raShis = [stf(1).ray.rangeShifter];
                %sanity check range shifters
                if ~isscalar(unique([raShis.ID]))
                    matRad_cfg.dispError('MCsquare only supports one range shifter setting (on or off) per energy! Aborting.\n');
                end
                stfMCsquare(i).energyLayer(layerIx).rangeShifter = raShis(1);
            end

        end

        % remember order
        counterMCsquare = 0;
        MCsquareOrder = NaN * ones(dij.totalNumOfBixels,1);
        for i = 1:length(stf)
            for j = 1:size(MCsquareBDL.focusTable,1)
                for k = 1:numel(stfMCsquare(i).energyLayer(j).numOfPrimaries)
                    counterMCsquare = counterMCsquare + 1;
                    ix = find( ...
                        i                                         == dij.beamNum & ...
                        stfMCsquare(i).energyLayer(j).rayNum(k)   == dij.rayNum & ...
                        stfMCsquare(i).energyLayer(j).bixelNum(k) == dij.bixelNum ...
                        );

                    MCsquareOrder(ix) = counterMCsquare;
                end
            end
        end

        if any(isnan(MCsquareOrder))
            matRad_cfg.dispError('Invalid ordering of Beamlets for MCsquare computation!');
        end

        %% Write config files
        % override HU_Density_Conversion_File and HU_Material_Conversion_File in case of Heterogeneity density sampling
        if isfield(ct,'modulated') && ct.modulated
            % copy and override with default HU conversion files
            if ~isfolder(['Scanners' '/' 'densitySampling'])
                mkdir(['Scanners' '/' 'densitySampling'])
            end
            copyfile(pln.propMC.HU_Density_Conversion_File,['Scanners' '/' 'densitySampling' '/' 'HU_Density_Conversion.txt'])
            copyfile(pln.propMC.HU_Material_Conversion_File,['Scanners' '/' 'densitySampling' '/' 'HU_Material_Conversion.txt'])

            % prepare sampled densities and combine with HU
            sampledDensities(1,:) = 6000:5999+length(ct.sampledDensities);
            sampledDensities(2,:) = ct.sampledDensities;

            % write sampled densities
            fID = fopen(['Scanners' '/' 'densitySampling' '/' 'HU_Density_Conversion.txt'],'a');
            fprintf(fID,'\n%i       %.3f',sampledDensities);
            fclose(fID);

            % write material conversion
            fID = fopen(['Scanners' '/' 'densitySampling' '/' 'HU_Material_Conversion.txt'],'a');
            % fprintf(fID,'\n6000    40      # Schneider_Lung');
            %             fprintf(fID,'\n6000    17      # Water');
            fprintf(fID,'\n6000    14      # Lung (use 40 for Schneider_Lung)');

            fclose(fID);

            % set custom HU conversion files to be used by MCsquare
            pln.propMC.HU_Density_Conversion_File = ['Scanners' '/' 'densitySampling' '/' 'HU_Density_Conversion.txt'];
            pln.propMC.HU_Material_Conversion_File = ['Scanners' '/' 'densitySampling' '/' 'HU_Material_Conversion.txt'];
        end

        % write patient data
        MCsquareBinCubeResolution = [dij.doseGrid.resolution.x ...
            dij.doseGrid.resolution.y ...
            dij.doseGrid.resolution.z];
        pln.propMC.writeMhd(HUcube{ctScen},MCsquareBinCubeResolution);

        % Copy necessary files to MCrun folder
        copyfile(pln.propMC.HU_Density_Conversion_File, [pln.propMC.MCrun_Directory 'HU_Density_Conversion.txt'])
        pln.propMC.HU_Density_Conversion_File = [pln.propMC.MCrun_Directory 'HU_Density_Conversion.txt'];
        copyfile(pln.propMC.HU_Material_Conversion_File, [pln.propMC.MCrun_Directory 'HU_Material_Conversion.txt'])
        pln.propMC.HU_Material_Conversion_File = [pln.propMC.MCrun_Directory 'HU_Material_Conversion.txt'];
        copyfile(pln.propMC.BDL_Machine_Parameter_File, [pln.propMC.MCrun_Directory 'BDL_file.txt'])
        pln.propMC.BDL_Machine_Parameter_File = [pln.propMC.MCrun_Directory 'BDL_file.txt'];

        % write config file
        pln.propMC.writeMCsquareinputAllFiles(MCsquareConfigFile,stfMCsquare,MCsquareBDL.focusTable);

        if isfield(dij,'RBE_model')
            [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,ct.numOfCtScen);
        end

        % write parameters to a MCparam file that can be used to later read the dose back in
        MCparam.dij = dij; % this can be done here since the dij is not filled at this point
        if isfield(dij,'physicalDose')
            MCparam.dij = rmfield(MCparam.dij,'physicalDose');
        end
        if isfield(dij,'mLETDose')
            MCparam.dij = rmfield(MCparam.dij,'mLETDose');
        end
        MCparam.VdoseGrid = VdoseGrid;
        MCparam.calcDoseDirect = calcDoseDirect;
        MCparam.totalWeights = sum(totalWeights);
        MCparam.Beamlet_Mode = pln.propMC.Beamlet_Mode;
        MCparam.nbHistoriesTotal = pln.propMC.numHistories;
        MCparam.MCsquareOrder = MCsquareOrder;

        % Generate output folder and save MCparam
        if ~exist(pln.propMC.Output_Directory,'dir')
            mkdir(pln.propMC.Output_Directory)
        end
        save([pln.propMC.Output_Directory '/' 'MCparam.mat'],'MCparam')

        %% MC computation and dij filling
        % run MCsquare
        if ~pln.propMC.externalCalculation
            mcSquareCall = [mcSquareBinary ' ' MCsquareConfigFile];
            matRad_cfg.dispInfo(['Calling Monte Carlo Engine: ' mcSquareCall]);
            [status,cmdout] = system(mcSquareCall,'-echo');
        end

        % Skip readout if external files were generated
        if ~pln.propMC.externalCalculation
            dij = pln.propMC.readFiles(strrep(pln.propMC.Output_Directory, '/', filesep));
        else
            dij = struct([]);
        end

        if ~pln.propMC.externalCalculation
            matRad_cfg.dispInfo('Simulation finished!\n');
        else
            matRad_cfg.dispInfo('Files generated for external calculation!\n');
        end

    end

    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(scenarioIx,:);
    end
end

% Order fields for easier comparison between different dijs
dij = orderfields(dij);

%% cd back
cd(currFolder);

if ishandle(figureWait)
    delete(figureWait);
end

end
