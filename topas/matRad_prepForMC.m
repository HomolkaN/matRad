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

cstHetero = matRad_cstHeteroAutoassign(cst);
carbHetero = matRad_calcDoseDirect(ct,stf,pln,cstHetero,resultGUI.w,param);
matRad_compareDose(resultGUI.physicalDose,carbHetero.physicalDose, ct, cst,[1 0 0]);
save('matfiles/Boxphantom_protons/heterogeneity/boxphantom_protons_homogeneous.mat','resultGUI','carbHetero','pln','cst','ct','stf')


%%
% sigma = 0.24;
% variance = sigma^2;
% 
% upper = 1;
% lower = 0;
% mu = 0.4;
% 
% pd = makedist('beta');
% pd.a = mu    * ( mu*(1-mu)/variance - 1 );
% pd.b = (1-mu)* ( mu*(1-mu)/variance - 1 );
% %%
% num = length(ct.cube{1}(cst{3,4}{1}));
% 
% for n = 1:100
%     clear X gwn
%     X = zeros(num,1);
%     for i = 1:num
%         X(i) = pd.random;
%     end
%     ct.cube{1}(cst{3,4}{1}) = X;
%     ct.cubeHU{1}(cst{3,4}{1}) = -1000 + 1000 * X;
%     save(['matfiles/Boxphantom_protons/heterogeneity/files/data_' num2str(n,'%03.f') '.mat'],'resultGUI','pln','cst','ct','stf')
% end


