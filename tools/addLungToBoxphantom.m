clear cst test*
load BOXPHANTOM

%%
cst{3,1} = 2;
cst{3,2} = 'LungBox';
cst{3,3} = 'OAR';
%%
cst{3,5}.TissueClass = 1;
cst{3,5}.alphaX = 0.1;
cst{3,5}.betaX = 0.05;
cst{3,5}.Priority = 2;
cst{3,5}.Visible = true;
cst{3,5}.visibleColor = [0 0.3333 0.6667];
%%
cst{3,6}.type = 'square overdosing';
cst{3,6}.dose = 42;
cst{3,6}.penalty = 400;
cst{3,6}.EUD = NaN;
cst{3,6}.volume = NaN; 
cst{3,6}.robustness = 'none';


%% Setting the coordinates of the box
% centerP = [90,240,240];
% height = 60;
% width = 240;
% depth = 240;
% %%
% resolution(1)= ct.resolution.x;
% resolution(2)= ct.resolution.y;
% resolution(3)= ct.resolution.z;
% 
% centerP_corr = centerP ./ resolution;
% height_corr = height/resolution(1);
% width_corr = width/resolution(2);
% depth_corr = depth/resolution(3);
centerP_corr = [30,80,80];
height_corr = 20;
width_corr = 80;
depth_corr = 80;


mask = zeros(ct.cubeDim);

for i=-height_corr/2:height_corr/2
    for j=-width_corr/2:width_corr/2
        for k=-depth_corr/2:depth_corr/2
            mask(centerP_corr(1)+i,centerP_corr(2)+j,centerP_corr(3)+k) = 1;
        end
    end
end
cst{3,4}{1} = find(mask == 1);



%% move target 10mm
clear mask center* height* width* depth*
centerP_corr = [70,80,80];
height_corr = 20;
width_corr = 20;
depth_corr = 20;


mask = zeros(ct.cubeDim);

for i=-height_corr/2:height_corr/2
    for j=-width_corr/2:width_corr/2
        for k=-depth_corr/2:depth_corr/2
            mask(centerP_corr(1)+i,centerP_corr(2)+j,centerP_corr(3)+k) = 1;
        end
    end
end
cst{2,4}{1} = find(mask == 1);


%% Creating the cube in the CT
vIxLung = cst{3,4}{1};
ct.cubeHU{1}(vIxLung) = -600;
%ct.cubeHU{1}(cst{3,4}{1,1}) = 1;
clearvars -except ct cst

%%
save BOXPHANTOM_LUNG_CORR

%%
% 
figure
plane = 3;
slice = 80;
matRad_plotCtSlice(gca,ct.cubeHU,1,plane,slice)
matRad_plotVoiContourSlice(gca,cst,ct.cubeHU,1,[],plane,slice)
