function matRad_heterogeneityCorrectionWrapper(patient,particle,machine,doseMode,param)
% matRad wrapper for heterogeneity testing with available preimported
% patients
%
% call
%   matRad_heterogeneityCorrectionWrapper(patient,mode,machine,doseMode,param)
%
% input
%   patient:         preloaded patient as string
%                    ('S001','S006','H3368','Boxphantom')
%   mode:            particle mode ('protons' or 'carbon')
%   machine:         machine file as string ('generic' or 'HIT_APM')
%   doseMode:        dose mode (physicalDose or RBE) as string
%   param (optional)
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
%%
if exist('param','var')
    if ~isfield(param,'logLevel')
        param.logLevel = 1;
    end
else
    param.subIx          = [];
    param.logLevel       = 1;
end

%% Patient Data Import
% if contains(patient,'006')
%     load('dataImport/lowRes/S00006-ID20180322.mat');
% elseif contains(patient,'001')
%     load('dataImport/lowRes/S00001_ID-20171201.mat');
% elseif contains(patient,'3368')
%     load('dataImport/lowRes/H03368_ID-20171201.mat');
% else
load(['C:\Home\homolka\git\matRad\dataImport\',patient,'.mat'])

if contains(patient,'Boxphantom')
    load('BOXPHANTOM_LUNG_PLN.mat');
    cst{2,6}.dose = 30;
    pln.numOfFractions = 6;
end

cst = matRad_slimCST(cst);
%%
pln.radiationMode = particle;
pln.machine       = machine;

if contains(doseMode,'physical')
    modelName    = 'none';
    quantityOpt  = 'physicalDose';
    
elseif contains(doseMode,'RBE')
    quantityOpt  = 'RBExD';
    
    if contains(particle,'carbon')
        modelName    = 'LEM';
    elseif contains(particle,'helium')
        modelName    = 'HEL';
    elseif contains(particle,'proton')
        if contains(doseMode,'MCN')
            modelName    = 'MCN';
        elseif contains(doseMode,'const')
            modelName = 'constRBE';
        end
    end
else
    error('Must choose valid dose cube');
end


%%
% The remaining plan parameters are set like in the previous example files
%pln.numOfFractions        = 30;
%pln.propStf.gantryAngles  = 0;
%pln.propStf.couchAngles   = 0;
if isnan(pln.propStf.bixelWidth) || isempty(pln.propStf.bixelWidth)
    pln.propStf.bixelWidth    = 3;
end
%pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
%pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propOpt.runDAO        = 0;
%pln.propOpt.runSequencing = 0;

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
pln.heterogeneity.useDoseCurves = false;
pln.heterogeneity.useOrgDepths = false;


%% Dose Calculation
% The previously set parameters will be automatically transferred into the
% dose calculation and checked for consistency. If nothing is specified, it
% will automatically be disabled.
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
clear dij
%% Visualize differences
% try
%     if contains(doseMode,'Physical')
%         [exp1,exp2,fig] = matRad_compareDose(carbHomo.physicalDose,carbHetero.physicalDose,ct,cst,[1 0 0]);
%     else
%         [exp1,exp2,fig] = matRad_compareDose(carbHomo.RBExD,carbHetero.RBExD,ct,cst,[1 0 0]);
%     end
%     savefig(fig.axial.fig,[patient,'_',mode,'_',machine,'_',doseMode,'.fig'])
% catch
%     warning(['Could not compare dose cubes']);
% end

save(['Results\',patient,'_',modelName,'_',particle,'_',machine,'_',doseMode,'.mat'])

end