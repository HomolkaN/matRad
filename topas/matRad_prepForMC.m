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
matRad_rc

%% Patient Data Import
%load('BOXPHANTOM_LUNG_NARROW_save');
load('BOXPHANTOM_LUNG_NARROW_1mm');

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
pln.radiationMode = 'protons';
pln.machine       = 'generic_TOPAS_cropped_APM';
% pln.machine       = 'GenericTest';

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

modelName           = 'none';
quantityOpt         = 'physicalDose';

%modelName           = 'none';
%quantityOpt         = 'PhysicalDose';

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions = 6;
%pln.propStf.lateralSpotSpacing = 40;
%pln.propStf.longitudinalSpotSpacing = 50;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 2;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
% pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
% pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

pln.propDoseCalc.doseGrid.resolution.x = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 1; % [mm]


pln.propDoseCalc.resolution = ct.resolution;


% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario

%% Generate Beam Geometry STF
stf1 = matRad_generateStf(ct,cst,pln,param);
stf = stf1;

%% Selecting the center ray
for r = 1:length(stf1.ray)
    if stf1.ray(r).rayPos_bev(1) == 0 && stf1.ray(r).rayPos_bev(2) == 0 && stf1.ray(r).rayPos_bev(3) == 0
        break
    end
end
stf.totalNumOfBixels = 1;
stf.ray = stf1.ray(r);
energyIx = round(numel(stf.ray.energy)/2);
stf.ray.energy = stf.ray.energy(energyIx);
stf.ray.focusIx = stf.ray.focusIx(energyIx);
stf.ray.rangeShifter = stf.ray.rangeShifter(energyIx);
stf.numOfRays = 1;
stf.numOfBixelsPerRay = 1;


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

%% Dose Calculation
% The previously set parameters will be automatically transferred into the
% dose calculation and checked for consistency. If nothing is specified, it
% will automatically be disabled.
dij = matRad_calcParticleDose(ct,stf,pln,cst,param);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
resultGUI = matRad_fluenceOptimization(dij,cst,pln,param);
%resultGUI = rmfield(resultGUI,'wUnsequenced');
%resultGUI = rmfield(resultGUI,'physicalDose_beam1');
%resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
%%
% save('matfiles/Boxphantom_protons/heterogeneity/boxphantom_protons1mm_homogeneous.mat','resultGUI','pln','cst','ct','stf')
cstHetero = matRad_cstHeteroAutoassign(cst);
resultGUIhetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w,param);
%matRad_compareDose(resultGUI.physicalDose,hetero_physicalDose.physicalDose, ct, cst,[1 0 0]);
save('boxphantom_protons1mm_inhomogeneous.mat','resultGUIhetero','pln','cst','ct','stf')

%%
p = 0.27;
variance = p*(1-p);
display(['Variance is var = ' num2str(variance)]);
sigma = sqrt(variance);

upper = 1;
lower = 0;
mu = 0.4;

% dist = makedist('beta');
% dist.a = mu    * ( mu*(1-mu)/variance - 1 );
% dist.b = (1-mu)* ( mu*(1-mu)/variance - 1 );

pd = makedist('Normal');
pd.mu = mu;
pd.sigma = sigma;
dist = truncate(pd,lower,upper);

%
X = -0.1:0.001:1.1;
f = figure;
pDist = pdf(dist,X);
plot(X,pDist/max(pDist))
xlim([-0.02 1.02])
ylim([0 1.1])
xlabel('Dichte des Voxels')
ylabel('PDF')
title(['var = ' num2str(variance) ', p = ' num2str(p)])

%%
% X = 0:0.01:1;
% variance = X.*(1-X);
% figure
% plot(X,variance)

%%
num = length(ct.cube{1}(cst{3,4}{1}));

for n = 1:100
    clear X gwn
    X = random(dist,1,num);
    
    %Bernoulli
    %X = rand(num,1);
    %X = double(X < mu);
    
    ct.cube{1}(cst{3,4}{1}) = X;
    ct.cubeHU{1}(cst{3,4}{1}) = -1000 + 1000 * X;
%   save(['matfiles/Boxphantom_protons/heterogeneity/files/data_' num2str(n,'%03.f') '.mat'],'resultGUI','pln','cst','ct','stf')
    save(['data/data_' num2str(n,'%03.f') '.mat'],'resultGUI','pln','cst','ct','stf')
end


