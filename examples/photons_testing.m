%% Example: Photon Treatment Plan using VMC++ dose calculation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation based on the VMC++ Monte Carlo algorithm 
% (iii) how to inversely optimize the beamlet intensities directly from command window in MATLAB. 
% (iv) how to visualize the result

%% set matRad runtime configuration
matRad_rc %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the boxphantom
% into your workspace. 
load('BOXPHANTOM.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.

pln.radiationMode           = 'photons';  
pln.machine                 = 'Generic';
pln.numOfFractions          = 30;
pln.propStf.gantryAngles    = 0;
pln.propStf.couchAngles     = 0;
pln.propStf.bixelWidth      = 10;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runSequencing   = 0;
pln.propOpt.runDAO          = 0;

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

% retrieve bio model parameters
pln.bioParam = matRad_BioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
% stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities using 
% a Monte Carlo algorithm
% dij = matRad_calcPhotonDose(ct,stf,pln,cst);
% resultGUI = matRad_fluenceOptimization(dij,cst,pln);
% 
% [dij2,stf2] = matRad_spotRemoval(dij,resultGUI.w,stf,0.03);
% resultGUI = matRad_fluenceOptimization(dij2,cst,pln);
% 
% %% Inverse Optimization for IMRT
% pln.propMC.numHistories = 1e4;
% 
% tic
% pln.propMC.engine = 'ompMC';
% resultGUI_ompMC = matRad_calcDoseDirectMC(ct,stf2,pln,cst,resultGUI.w);
% time_ompMC = toc;
% %%
% 
% tic
% pln.propMC.engine = 'TOPAS';
% % pln.propMC.beamProfile = 'uniform';
% pln.propMC.beamProfile = 'phasespace';
% % pln.propMC.beamProfile = 'virtualGaussian';
% resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf2,pln,cst,resultGUI.w);
% time_topas = toc;
% %%
% matRad_compareDose(resultGUI.physicalDose,resultGUI_ompMC.physicalDose,ct,cst,[1 1 1],1);

%%
stf = matRad_generateStfSinglePencilBeam(ct,cst,pln);
stf2 = matRad_generateStf(ct,cst,pln);
stf.ray = stf2.ray(find(sum(ismember(vertcat(stf2.ray.rayPos),[0 0 0]),2)==3));

resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,ones(stf.totalNumOfBixels,1));

%%
pln.propMC.numHistories = 1e4;

tic
pln.propMC.engine = 'ompMC';
resultGUI_ompMC = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w);
time_ompMC = toc;

tic
pln.propMC.engine = 'TOPAS';
pln.propMC.beamProfile = 'uniform';
resultGUI_TOPAS = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w);
time_topas = toc;


matRad_compareDose(resultGUI.physicalDose,resultGUI_ompMC.physicalDose,ct,cst,[1 1 1],1);



