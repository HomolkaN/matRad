function [pln,cst] = matRad_setPlanParameters(ct,cst,baseData,RBE_model,testing,pln)
if ~exist('pln','var')
    pln = struct;
end
if nargin < 5
    testing = false;
end
if ~isfield(pln,'propStf') || ~isfield(pln.propStf,'couchAngles')
    pln.propStf.gantryAngles = 0;
    pln.propStf.couchAngles = 0;
end
if testing

    pln.propStf.bixelWidth = 5;
else
    pln.propStf.bixelWidth = 3;
end
pln.propStf.numOfBeams = size(pln.propStf.gantryAngles,2);
pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*matRad_getIsoCenter(cst,ct);
pln.propStf.longitudinalSpotSpacing = 1;

pln.numOfFractions = 30;

target = find(cellfun(@(teststr) ~isempty(strfind(lower(teststr),lower('TARGET'))), cst(:,3)));
if ~iscell(cst{target(1),6})
    cst = matRad_computeVoiContoursWrapper(cst,ct);
end
if ~iscell(cst{target(1),6})
    cst = matRad_convertOldCstToNewCstObjectives(cst);
end
for i = target'
    cst{i,6}{1}.parameters{1} = 60;
end
if size(cst,2)>6 && iscell(cst{1,7})
    cst(:,7)=[];
end

s = split(baseData,'_');
pln.radiationMode   = s{1};
pln.machine = strjoin(s(2:end),'_');

if any(cellfun(@(teststr) ~isempty(strfind(lower(RBE_model),lower(teststr))), {'WED','MCN','LEM','RBE','HEL'}))
    quantityOpt                 = 'RBExD';
    pln.propOpt.bioOptimization = 'RBExD';
else
    quantityOpt                 = 'physicalDose';
    pln.propOpt.bioOptimization = 'physicalDose';
end

pln.propDoseCalc.lateralCutOff = 0.99;
pln.propDoseCalc.airOffsetCorrection = 1;

% retrieve scenarios for dose calculation and optimziation
% pln = rmfield(pln,{'bioParam','multScen','DicomInfo'});
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,RBE_model);
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario
if testing
    pln.propDoseCalc.resolution = struct('x',3,'y',3,'z',3);
else
    pln.propDoseCalc.resolution = struct('x',3,'y',3,'z',3);
end

pln.propHeterogeneity.calcHetero = false;
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% Monte Carlo settings
if testing
    pln.propMC.numHistories = 1e4;
    pln.propMC.numOfRuns = 1;
    pln.propMC.externalCalculation = false;
else
    pln.propMC.numHistories = 1e8;
    pln.propMC.numOfRuns = 5;
    pln.propMC.externalCalculation = true;
end
pln.propMC.engine = 'TOPAS';

end