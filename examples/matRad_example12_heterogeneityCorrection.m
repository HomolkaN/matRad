%% Example: Heterogeneity corrected carbon Ion Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a carbon ion dose calculation plan including variable RBE optimization
% (iii) how to inversely optimize the pencil beam intensities based on the
% RBE-weighted dose
% (iv) how to inversely optimize the pencil beam intensities based on the
% biological effect
% (v) how to change the tissues' radiobiological characteristics
% (vi) how to recalculated the dose considering the previously optimized pencil beam intensities
% (vii) how to compare the two results

%% set matRad runtime configuration
clear

param.calcDoseDirect = false;
param.subIx = [];
param.logLevel = 1;

%% Patient Data Import
load('BOXPHANTOM_LUNG');
cstDev = cst;

for i = 1:3
cst{i,6} = struct;
cst{i,6}.dose  = cstDev{i,6}{1,1}.parameters{1,1};
cst{i,6}.penalty = cstDev{i,6}{1,1}.penalty;
cst{i,6}.robustness = 'none';
end

cst{1,6}.type = 'square overdosing';
cst{2,6}.type = 'square deviation';  
cst{3,6}.type = 'square overdosing';
clear cstDev
%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use carbon ions for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. In order to use heterogeneity correction, the base data must contain
% the depth information as a struct. For this purpose, matRad features generic base data in the file
% 'carbon_GenericAPM.mat'; consequently the machine has to be set accordingly

pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'generic_TOPAS_cropped_APM';

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
%('none': physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN': McNamara-variable RBE model for protons; 
% 'WED':  Wedenberg-variable RBE model for protons
% 'LEM': local effect model 
% As we use carbons, we use the local effect model.
% Therefore we set modelName to LEM

modelName           = 'LEM';
quantityOpt         = 'RBExD';   

%modelName           = 'none';
%quantityOpt         = 'PhysicalDose';   

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions = 6;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario                                            

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln,param);

%% Setting heterogeneity parameters
% To enable heterogeneity correction, the parameter calcHetero must be set
% to true. The type of correction ('complete', 'depthBased' or 'voxelwise')
% can be set. It will be set to 'complete' by default, if nothing is
% specified. Independently, the parameter useDoseCurves enables the
% calculation of RBE using fitted alpha and sqrtBeta curves implemented in
% the APM base data files. 

pln.heterogeneity.calcHetero = true;
pln.heterogeneity.type = 'complete'; % optional
pln.heterogeneity.useDoseCurves = true;
pln.heterogeneity.useOrgDepths = false;

%%

stf1 = stf;

%% Selecting the center ray
for r = 1:length(stf1.ray)
    if stf1.ray(r).rayPos_bev(1) == 0 && stf1.ray(r).rayPos_bev(2) == 0 && stf1.ray(r).rayPos_bev(3) == 0
        break
    end
end

% for t = 1:length(stf1.ray)
%     if stf1.ray(t).rayPos_bev(1) == -10 && stf1.ray(t).rayPos_bev(2) == 0 && stf1.ray(t).rayPos_bev(3) == 0
%         break
%     end
% end


stf.ray = stf1.ray(r);
energyIx = round(numel(stf1.ray(r).energy)/2);
%stf.ray.energy = stf1.ray(1).energy(energyIx);
stf.ray.energy = stf1.ray(r).energy(energyIx);
stf.ray.focusIx = stf1.ray(r).focusIx(energyIx);
stf.ray.rangeShifter = stf1.ray(r).rangeShifter(energyIx);

% stf.ray(2) = stf1.ray(t);
% energyIx = round(numel(stf1.ray(t).energy)/2);
% stf.ray(2).energy = stf1.ray(t).energy(energyIx);
% stf.ray(2).focusIx = stf1.ray(t).focusIx(energyIx);
% stf.ray(2).rangeShifter = stf1.ray(t).rangeShifter(energyIx);
stf.numOfRays = 1;
stf.numOfBixelsPerRay = 1;
stf.totalNumOfBixels = 1;

%%
dij = matRad_calcParticleDose(ct,stf,pln,cst,param);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
carbHomo = matRad_fluenceOptimization(dij,cst,pln,param);

%% Assign heterogeneity flags for lung tissue
% The above calculation will give you a warning, specifying that the cst
% file does not contain heterogeneity information. This information is
% needed in order to correctly calculate through which organs the
% correction has to be considered. In order to change this and assign
% heterogeneity correction to the cst file, we can use the
% cstHeteroAutoassign function. This will automatically specify lung tissue
% with the heterogeneity flag for 'lung'.
cstHetero = matRad_cstHeteroAutoassign(cst);

%% Calculate dose again with heterogeneityCorrection
carbHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,carbHomo.w,param);

%% Visualize differences
% matRad_compareDose(carbHomo.physicalDose,carbHetero.physicalDose,ct,cst,[1 0 0]);
matRad_compareDose(carbHomo.RBExD,carbHetero.RBExD,ct,ct,cst,[1 1 0]);


