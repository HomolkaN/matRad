function bixel = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData, heteroCorrDepths, propHeterogeneity , vTissueIndex)
% matRad visualization of two-dimensional dose distributions
% on ct including segmentation
%
% call
%   dose = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
%   baseData:       base data required for particle dose calculation
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% skip heterogeneity correction for other functions
if nargin < 5
    heteroCorrDepths = [];
    % Load heterogeneity config for Gauss functions that are called even if heterogeneity correction is turned off
    propHeterogeneity = matRad_HeterogeneityConfig();
end

% Check if correct base data is loaded for heterogeneity correction
if ~isempty(heteroCorrDepths) && ~isstruct(baseData.Z) && ~strcmp(propHeterogeneity.type,'numerical')
    matRad_cfg.dispWarning('calcParticleDoseBixel: heterogeneity correction enabled but no APM base data was loaded.')
end

% add potential offset
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

%% interpolate depth dose, sigmas and weights and calculate lateral sigmas
% This is not heterogeneity corrected
if isfield(baseData,'sigma1')
    
    % interpolate depth dose, sigmas, and weights
    X = matRad_interp1(depths,[baseData.Z baseData.weight baseData.sigma1 baseData.sigma2],radDepths,'extrap');
    
    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,3).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;
    
else
    
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,[baseData.Z baseData.sigma],radDepths);
    
    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;
    
    %compute lateral sigma
    sigmaSq = X(:,2).^2 + sigmaIni_sq;
    
end

%  calculate lateral profiles
if isfield(baseData,'sigma1')
    
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr); % Gauss
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    
    bixel.L = baseData.LatCutOff.CompFac * ((1-X(:,2)).*L_Narr + X(:,2).*L_Bro); % (1-w)*L_Narr + w*L_Bro
    
else
    
    bixel.L = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) ./(2*pi*sigmaSq);
    
end

%% calculate sigma in range direction
bixel.Z = X(:,1);

%% calculating the physical dose
bixel.physDose = conversionFactor * bixel.L .* bixel.Z;

%% check if we have valid dose values
if any(isnan(bixel.physDose)) || any(bixel.physDose<0)
    error('Error in particle dose calculation.');
end
