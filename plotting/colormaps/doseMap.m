function colorMap = doseMap(cMapSize)
% matRad difference colormap
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


colorMapData = ...
    [0	    0	    0.6000
    0	    0	    0.7000
    0	    0	    0.8000
    0	    0	    0.9000
    0	    0	    1
    0	    0.1000	1
    0	    0.2000	1
    0	    0.3000	1
    0	    0.4000	1
    0	    0.5000	1
    0	    0.6000	1
    0	    0.7000	1
    0	    0.8000	1
    0	    0.9000	1
    0	    1	    1
    0.1000	1	    0.9000
    0.2000	1	    0.8000
    0.3000	1	    0.7000
    0.4000	1	    0.6000
    0.5000	1	    0.5000
    0.6000	1	    0.4000
    0.7000	1	    0.3000
    0.8000	1	    0.2000
    0.9000	1	    0.1000
    1	    1	    0
    1	    0.9000	0
    1	    0.8000	0
    1	    0.700000000000000	0
    1	    0.600000000000000	0
    1	    0.500000000000000	0
    1	    0.400000000000000	0
    1	    0.300000000000000	0
    1	    0.200000000000000	0
    1	    0.100000000000000	0
    1	    0	0
    0.900000000000000	0	0
    0.800000000000000	0	0
    0.700000000000000	0	0
    0.600000000000000	0	0
    0.500000000000000	0	0
    0.500000000000000	0	0
    0.428571428571429	0	0
    0.357142857142857	0	0
    0.285714285714286	0	0
    0.214285714285714	0	0
    0.142857142857143	0	0
    0.0714285714285715	0	0
    0	0	0
    0	0	0
    0.0833333333333333	0.0520800000000000	0.0331666666666667
    0.166666666666667	0.104160000000000	0.0663333333333333
    0.250000000000000	0.156240000000000	0.0995000000000000
    0.333333333333333	0.208320000000000	0.132666666666667
    0.416666666666667	0.260400000000000	0.165833333333333
    0.500000000000000	0.312480000000000	0.199000000000000
    0.583333333333333	0.364560000000000	0.232166666666667
    0.666666666666667	0.416640000000000	0.265333333333333
    0.750000000000000	0.468720000000000	0.298500000000000
    0.833333333333333	0.520800000000000	0.331666666666667
    0.916666666666667	0.572880000000000	0.364833333333333
    1	0.624960000000000	0.398000000000000
    1	0.677040000000000	0.431166666666667
    1	0.729120000000000	0.464333333333333
    1	0.781200000000000	0.497500000000000];

if nargin < 1
    colorMap = colorMapData;
elseif size(colorMapData,1) == cMapSize
    colorMap = colorMapData;
else
    %We have to interpolate the colormap
    newX = linspace(1,64,cMapSize);
    oldX = 1:64;
    colorMap = interp1(oldX,colorMapData,newX);
    %{
    %resample via HSV.. more color-true than above, but doesn't work with
    %every colormap
    hsv                        = rgb2hsv(cm);
    hsv(144:end,1)             = hsv(144:end,1)+1;
    ColorMap                   = interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,cMapSize));
    ColorMap(cm(:,1)>1,1) = ColorMap(cm(:,1)>1,1)-1;
    ColorMap                   = hsv2rgb(cm);
    %}
end

end

