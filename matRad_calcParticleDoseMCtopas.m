function dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect,exportForExternalCalculation)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad TOPAS Monte Carlo proton dose calculation wrapper
%   This calls a TOPAS installation (not included in matRad due to
%   licensing model of TOPAS) for MC simulation
%
% call
%   dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel               number of histories per beamlet (nCasePerBixel > 1),
%                               max stat uncertainity (0 < nCasePerBixel < 1)
%   calcDoseDirect:             binary switch to enable forward dose
%                               calcualtion
% output
%   dij:                        matRad dij struct
%
% References
%
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

matRad_cfg = MatRad_Config.instance();

if nargin < 5
    % set number of particles simulated per pencil beam
    nCasePerBixel = matRad_cfg.propMC.particles_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end

if nargin < 6
    calcDoseDirect = false;
end

if nargin < 7
    exportForExternalCalculation = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    matRad_cfg.dispWarning('Variance scoring for TOPAS not yet supported.');
end

if isfield(pln,'propMC') && isfield(pln.propMC,'config')        
    if isa(pln.propMC.config,'MatRad_TopasConfig')
        matRad_cfg.dispInfo('Using given Topas Configuration in pln.propMC.config!\n');
        topasConfig = pln.propMC.config;
    else 
        %Create a default instance of the configuration
        topasConfig = MatRad_TopasConfig();
        
        %Overwrite parameters
        %mc = metaclass(topasConfig); %get metaclass information to check if we can overwrite properties
        
        if isstruct(pln.propMC.config)
            props = fieldnames(pln.propMC.config);
            for fIx = 1:numel(props)
                fName = props{fIx};
                if isprop(topasConfig,fName)
                    %We use a try catch block to catch errors when trying
                    %to overwrite protected/private properties instead of a
                    %metaclass approach
                    try 
                        topasConfig.(fName) = pln.propMC.config.(fName);
                    catch
                        matRad_cfg.dispWarning('Property ''%s'' for MatRad_TopasConfig will be omitted due to protected/private access or invalid value.',fName);
                    end
                else
                    matRad_cfg.dispWarning('Unkown property ''%s'' for MatRad_TopasConfig will be omitted.',fName);
                end
            end
        else
            matRad_cfg.dispError('Invalid Configuration in pln.propMC.config');
        end
    end
else
    topasConfig = MatRad_TopasConfig();
end
        

if ~calcDoseDirect
    matRad_cfg.dispWarning('You have selected TOPAS dij calculation, this may take a while ^^');
    pln.propMC.calcDij = true;
end

if ~isfield(pln.propStf,'useRangeShifter') 
    pln.propStf.useRangeShifter = false;
end

env = matRad_getEnvironment();

%% Initialize dose Grid as usual
% for TOPAS we explicitly downsample the ct to the dose grid (might not
% be necessary in future versions with separated grids)
matRad_calcDoseInit;
[ctR,~] = matRad_resampleTopasGrid(ct,cst,pln,stf);
% overwrite CT grid in dij in case of modulation.
if isfield(ctR,'ctGrid')
    dij.ctGrid = ctR.ctGrid;
end

% fill bixels, rays and beams in case of dij calculation
%if ~calcDoseDirect
counter = 1;
for f = 1:dij.numOfBeams
    for r = 1:stf(f).numOfRays
        for b = 1:stf(f).numOfBixelsPerRay(r)
            dij.bixelNum(counter) = b;
            dij.rayNum(counter)   = r;
            dij.beamNum(counter)  = f;
            counter = counter + 1;
        end
    end
end
%end

%% sending data to topas

load([pln.radiationMode,'_',pln.machine],'machine');
topasConfig = MatRad_TopasConfig();
% Create Base Data
topasConfig.radiationMode = stf.radiationMode;

topasBaseData = MatRad_TopasBaseData(machine,stf);%,TopasConfig);

topasConfig.numHistories = nCasePerBixel;
if isfield(pln.propMC,'numOfRuns')
    topasConfig.numOfRuns = pln.propMC.numOfRuns;
end
if exportForExternalCalculation
    if isfield(pln,'patientID')
        topasConfig.workingDir = [topasConfig.workingDir pln.radiationMode filesep pln.patientID '_'];
    end
    topasConfig.workingDir = [topasConfig.workingDir pln.machine,'_',pln.radiationMode];
end
% topasConfig.numOfRuns = matRad_cfg.propMC.topas_defaultNumBatches;

%Collect weights
if calcDoseDirect
    w = zeros(sum([stf(:).totalNumOfBixels]),ctR.numOfCtScen);
    counter = 1;
    for i = 1:length(stf)
        for j = 1:stf(i).numOfRays
            rayBix = stf(i).numOfBixelsPerRay(j);
            w(counter:counter+rayBix-1,:) = stf(i).ray(j).weight;
            counter = counter + rayBix;
        end
    end
end

if isfield(pln,'bioParam') && strcmp(pln.bioParam.quantityOpt,'RBExD')
    topasConfig.scorer.RBE = true;
    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
    dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
end

currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end    
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                
                %Overwrite CT (TEMPORARY - we should use 4D calculation in
                %TOPAS here)
%                 ctR.cubeHU = cubeHUresampled(ctScen);
%                 ctR.cube = cubeResampled(ctScen);
                %Delete previous topas files
                files = dir([topasConfig.workingDir,'*']);
                files = {files(~[files.isdir]).name};
                fclose('all');
                for i = 1:length(files)
                    delete([topasConfig.workingDir,files{i}])
                end
               
                if calcDoseDirect
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData,w(:,ctScen));
                else
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData);
                end
                
                % Run simulation for current scenario
                cd(topasConfig.workingDir);
                              
                if exportForExternalCalculation
                    save('dij.mat','dij')
                    save('weights.mat','w')
                    matRad_cfg.dispInfo('TOPAS simulation skipped for external calculation\n');
                else
                    
                    for beamIx = 1:numel(stf)
                        
                        for runIx = 1:topasConfig.numOfRuns
                            fname = sprintf('%s_field%d_run%d',topasConfig.label,beamIx,runIx);
                            if isfield(pln.propMC,'verbosity') && strcmp(pln.propMC.verbosity,'full')
                                topasCall = sprintf('%s %s.txt',topasConfig.topasExecCommand,fname);
                            else
                                topasCall = sprintf('%s %s.txt > %s.out > %s.log',topasConfig.topasExecCommand,fname,fname,fname);
                            end

                            if topasConfig.parallelRuns
                                finishedFiles{runIx} = sprintf('%s.finished',fname);
                                delete(finishedFiles{runIx});
                                topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
                            end
                            
                            matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
                            [status,cmdout] = system(topasCall,'-echo');
                            if status == 0
                                matRad_cfg.dispInfo('TOPAS simulation completed succesfully\n');
                            else
                                if status == 139
                                    matRad_cfg.dispError('TOPAS segmentation fault might be caused from an outdated TOPAS version or Linux distribution');
                                else
                                    matRad_cfg.dispError('TOPAS simulation exited with error code %d\n',status);
                                end
                            end
                        end
                        
                        if topasConfig.parallelRuns
                            runsFinished = false;
                            pause('on');
                            while ~runsFinished
                                pause(1);
                                fin = cellfun(@(f) exist(f,'file'),finishedFiles);
                                runsFinished = all(fin);
                            end
                        end
                        
                    end
                    
                    cd(currDir);
                    
                    %% Simulation finished - read out volume scorers from topas simulation
                    if calcDoseDirect
                        topasCubes = matRad_readTopasData(topasConfig.workingDir);
                    else
                        topasCubes = matRad_readTopasData(topasConfig.workingDir,dij);
                    end
                    
                    fnames = fieldnames(topasCubes);
                    dij.MC_tallies = fnames;
                    
                    if calcDoseDirect
                        if isfield(topasCubes,'physicalDose')
                            for d = 1:length(stf)
                                dij.physicalDose{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_beam',num2str(d)]),[],1);
                            end
                        end
                        if isfield(topasCubes,'doseToWater')
                            for d = 1:length(stf)
                                dij.doseToWater{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['doseToWater_beam',num2str(d)]),[],1);
                            end
                        end    
                        if isfield(topasCubes,'alpha_beam1')
                            for d = 1:length(stf)
                                dij.alpha{ctScen,1}(:,d)           = reshape(topasCubes.(['alpha_beam',num2str(d)]),[],1);
                                dij.beta{ctScen,1}(:,d)            = reshape(topasCubes.(['beta_beam',num2str(d)]),[],1);
                                
                                [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
                                dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
                                
                                dij.mAlphaDose{ctScen,1}(:,d)      = dij.physicalDose{ctScen,1}(:,d) .* dij.alpha{ctScen,1}(:,d);
                                dij.mSqrtBetaDose{ctScen,1}(:,d)   = sqrt(dij.physicalDose{ctScen,1}(:,d)) .* dij.beta{ctScen,1}(:,d);
                            end
                        end
                        if isfield(topasCubes,'physicalDose_std_beam1')
                            for d = 1:length(stf)
                                dij.physicalDose_std{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_std_beam',num2str(d)]),[],1);
                            end
                        end        
                        if isfield(topasCubes,'LET')
                            for d = 1:length(stf)
                                dij.LET{ctScen,1}(:,d)    = reshape(topasCubes.(['LET_beam',num2str(d)]),[],1);
                                dij.mLETDose{ctScen,1}(:,d) = dij.physicalDose{ctScen,1}(:,d) .*  dij.LET{ctScen,1}(:,d);
                            end
                        end   
                    else
                        for f = 1:numel(fnames)
                            for d = 1:stf(f).totalNumOfBixels
                                dij.physicalDose{1}(:,d) = reshape(topasCubes.(fnames{f}){d},[],1);
                            end
                        end
                    end
                end
            end
        end
    end
end
    
    % manipulate isocenter back
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
    end
end
