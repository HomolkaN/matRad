function distalFallOff = matRad_distalFallOff(ct,cube)

% Define fineGrid resolution for interpolation
fineGridResolution = 0.0001;

% calculate IDD
idd = matRad_calcIDD(cube,0,'y');

% calculate max dose and its index
[maxDose, maxIdx] = max(idd);

% interpolate dose after the peak to fine grid
fineGrid = maxIdx:fineGridResolution:maxIdx+10;
fineDose = interp1(maxIdx:maxIdx+10,idd(maxIdx:maxIdx+10),fineGrid);

% find indices of 80% and 20% of the max dose
[~,r80] = min(abs(fineDose-maxDose*0.8));
[~,r20] = min(abs(fineDose-maxDose*0.2));

% Calculate distal falloff in mm
distalFallOff = (r20-r80)*fineGridResolution*ct.resolution.y;

end

