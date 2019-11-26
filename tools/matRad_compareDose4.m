function [gammaCube,gammaPassRate,hfig] = matRad_compareDose4(cube1, cube2, cube3, cube4, ct, cst,enable , contours, pln, criteria, n,localglobal)
% Comparison of two dose cubes in terms of gamma index, absolute and visual difference
%
% call
%    compareDose = matRad_compareDose(cube1,cube2,ct,cst,criteria,n,localglobal)
%
% input
%   cube1:         dose cube 1 as an M x N x O array
%   cube2:         dose cube 2 as an M x N x O array
%   ct:            ct struct with ct cube
%   cst:           list of interesting volumes inside the patient as matRad
%                  struct (optional, does not calculate gamma and DVH)
%   enable         (optional) specify if all sections are evaluated
%                  boolean 3x1 array
%                  [1 0 0]: evaluate only basic plots
%                  [0 1 0]: evaluate only line profiles
%                  [0 0 1]: evaluate only DVH
%   contours       (optional) specify if contours are plotted,
%                  'on' or 'off'
%   pln            (optional) specify BioModel for DVH plot
%   criteria:      [1x2] vector (optional) specifying the distance to agreement
%                  criterion; first element is percentage difference,
%                  second element is distance [mm], default [3 3]
%   n:             number of interpolations (optional). there will be 2^n-1
%                  interpolation points. The maximum suggested value is 3.
%                  default n=0
%   localglobal:   parameter to choose between 'global' and 'local'
%                  normalization (optional)
%
%
% output
%
%   gammaCube:      result of gamma index calculation
%   gammaPassRate:  rate of voxels passing the specified gamma criterion
%                   evaluated for every structure listed in 'cst'.
%                   note that only voxels exceeding the dose threshold are
%                   considered.
%   hfig:           Figure handle struct for all 3 figures and subplots,
%                   indexed by plane names
%
% References gamma analysis:
%   [1]  http://www.ncbi.nlm.nih.gov/pubmed/9608475
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

%% check if cubes consistent
if ~isequal(size(cube1),size(cube2))
    error('dose cubes must be the same size\n');
end

if ~exist('localglobal','var')
    localglobal = 'global';
end
if ~exist('n','var')
    n = 0;
end
if ~exist('criteria','var')
    criteria = [3 3];
end
if ~exist('cst','var') || isempty(cst)
    cst = [];
    skip = 1;
    fprintf('[\bWarning: cst not specified, skipping DVH, isocenter in the center of the doseCube]\b\n');
else
    skip = 0;
end
if ~exist('pln','var') || isempty(pln)
    pln = [];
end
if ~exist('contours','var') || isempty(contours)
    contours = true;
end
if exist('contours','var') || ~isempty(contours)
    if strcmp(contours,'on')
        contours = true;
    else
        contours = false;
    end
end

if (~exist('enable','var') || isempty(enable)) || strcmp(enable,'all')
    enable = [1 1 1];
end
if enable(1)==0
    gammaCube=[];
    gammaPassRate=[];
end

%% Calculate iso-center slices and resolution
if isempty(cst)
    s=size(cube1);
    isoCenter = [s(1)/2 s(2)/2 s(3)/2];
else
    isoCenter = matRad_getIsoCenter(cst,ct,0);
end

resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];

slicename = {round(isoCenter(2)./resolution(2)),round(isoCenter(1)./resolution(1)),round(isoCenter(3)./resolution(3))};
doseWindow = [0 max([cube1(:); cube2(:)])];
planename = {'coronal','sagittal','axial'};

%% Get the gamma cube
if enable(1)==1
    disp('Calculating gamma index cube...');
    if exist('criteria','var')
        relDoseThreshold = criteria(1); % in [%]
        dist2AgreeMm     = criteria(2); % in [mm]
    else
        dist2AgreeMm     = 3; % in [mm]
        relDoseThreshold = 3; % in [%]
    end
    
    [gammaCube,gammaPassRate] = matRad_gammaIndex(cube1,cube2,resolution,criteria,[],n,localglobal,cst);
    
    
    %%% Calculate absolute difference cube and dose windows for plots
    differenceCube1  = cube2-cube1;
    differenceCube2  = cube4-cube2;
    differenceCube3  = cube4-cube3;
    maxim = max([max(differenceCube1(cst{1,4}{1})),max(differenceCube2(cst{1,4}{1})),max(differenceCube3(cst{1,4}{1}))]);
    % maxim = max([max(differenceCube1(:)),max(differenceCube2(:)),max(differenceCube3(:))]);
    doseDiffWindow  = [-maxim +maxim];
    doseGammaWindow = [0 max(gammaCube(:))];
    relativeDifferenceCube = ( differenceCube1 ./ cube1 )*100;
    relativeDifferenceCube = relativeDifferenceCube(~isnan(relativeDifferenceCube));
    relativeDifference = round(relativeDifferenceCube(find(relativeDifferenceCube == max(abs(relativeDifferenceCube)))),2);
    %  relativeDifference = max(abs(differenceCube(:)))/max(max(abs(cube1(:))),max(abs(cube2(:))))*100;
    
    %% Load colormap for difference
    diffCMap = matRad_getColormap('diffPolar');
    
    %% Plot everything
    % Plot dose slices
    if contours == false
        cstHandle = [];
    else
        cstHandle = cst;
    end
    
    windowXY{1} = [10 70]/2;
    windowXY{2} = [31 130]/2;
    alpha = .8;
    for plane=3
        disp(['Plotting ',planename{plane},' plane...']);
        
        % Initialize Figure
        hfig.(planename{plane}).('fig') = figure('Renderer', 'painters', 'Position', [10 50 800 800]);
        set(gcf,'Color',[1 1 1]);
        
        % Plot Dose 1
        hfig.(planename{plane}).('cube1').Axes = subplot(2,4,1);
        [hfig.(planename{plane}).('cube1').CMap,...
            hfig.(planename{plane}).('cube1').Dose,...
            hfig.(planename{plane}).('cube1').Ct,...
            hfig.(planename{plane}).('cube1').Contour,...
            hfig.(planename{plane}).('cube1').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,cube1,plane,slicename{plane},windowXY,[],alpha,colorcube,jet,doseWindow,[],100);
        figtitle = get(gca,'title');
        figtitle = figtitle.String;
        
        % Plot Dose 2
        hfig.(planename{plane}).('cube2').Axes = subplot(2,4,2);
        [hfig.(planename{plane}).('cube2').CMap,...
            hfig.(planename{plane}).('cube2').Dose,...
            hfig.(planename{plane}).('cube2').Ct,...
            hfig.(planename{plane}).('cube2').Contour,...
            hfig.(planename{plane}).('cube2').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,cube2,plane,slicename{plane},windowXY,[],alpha,colorcube,jet,doseWindow,[],100);
        
        % Plot Dose 3
        hfig.(planename{plane}).('cube3').Axes = subplot(2,4,3);
        [hfig.(planename{plane}).('cube3').CMap,...
            hfig.(planename{plane}).('cube3').Dose,...
            hfig.(planename{plane}).('cube3').Ct,...
            hfig.(planename{plane}).('cube3').Contour,...
            hfig.(planename{plane}).('cube3').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,cube3,plane,slicename{plane},windowXY,[],alpha,colorcube,jet,doseWindow,[],100);
        
        % Plot Dose 4
        hfig.(planename{plane}).('cube4').Axes = subplot(2,4,4);
        [hfig.(planename{plane}).('cube4').CMap,...
            hfig.(planename{plane}).('cube4').Dose,...
            hfig.(planename{plane}).('cube4').Ct,...
            hfig.(planename{plane}).('cube4').Contour,...
            hfig.(planename{plane}).('cube4').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,cube4,plane,slicename{plane},windowXY,[],alpha,colorcube,jet,doseWindow,[],100);
        
        % Plot absolute difference
        hfig.(planename{plane}).('diff1').Axes = subplot(2,4,5);
        [hfig.(planename{plane}).('diff1').CMap,...
            hfig.(planename{plane}).('diff1').Dose,...
            hfig.(planename{plane}).('diff1').Ct,...
            hfig.(planename{plane}).('diff1').Contour,...
            hfig.(planename{plane}).('diff1').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,differenceCube1,plane,slicename{plane},windowXY,[],alpha,colorcube,diffCMap,doseDiffWindow,[],100);
        
        hfig.(planename{plane}).('diff2').Axes = subplot(2,4,7);
        [hfig.(planename{plane}).('diff2').CMap,...
            hfig.(planename{plane}).('diff2').Dose,...
            hfig.(planename{plane}).('diff2').Ct,...
            hfig.(planename{plane}).('diff2').Contour,...
            hfig.(planename{plane}).('diff2').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,differenceCube2,plane,slicename{plane},windowXY,[],alpha,colorcube,diffCMap,doseDiffWindow,[],100);
        
        hfig.(planename{plane}).('diff3').Axes = subplot(2,4,8);
        [hfig.(planename{plane}).('diff3').CMap,...
            hfig.(planename{plane}).('diff3').Dose,...
            hfig.(planename{plane}).('diff3').Ct,...
            hfig.(planename{plane}).('diff3').Contour,...
            hfig.(planename{plane}).('diff3').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cstHandle,1,differenceCube3,plane,slicename{plane},windowXY,[],alpha,colorcube,diffCMap,doseDiffWindow,[],100);
        
        set(hfig.(planename{plane}).('fig'),'name',figtitle);
        
        %% Adjusting axes
        
        % matRad_plotAxisLabels(hfig.(planename{plane}).('cube1').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube1').Axes, 'title'), 'string', 'matRad Homogen');
        
        % matRad_plotAxisLabels(hfig.(planename{plane}).('cube2').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube2').Axes, 'title'), 'string', 'matRad Heterogen');
        
        % matRad_plotAxisLabels(hfig.(planename{plane}).('cube3').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube3').Axes, 'title'), 'string', 'TOPAS Homogen');
        
        % matRad_plotAxisLabels(hfig.(planename{plane}).('cube4').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube4').Axes, 'title'), 'string', 'TOPAS Heterogen');
        
        
        %  matRad_plotAxisLabels(hfig.(planename{plane}).('diff1').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('diff1').Axes, 'title'), 'string', 'matRad Heterogen - matRad Homogen');
        
        %  matRad_plotAxisLabels(hfig.(planename{plane}).('diff2').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('diff2').Axes, 'title'), 'string', 'TOPAS - matRad Heterogen');
        
        %  matRad_plotAxisLabels(hfig.(planename{plane}).('diff3').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('diff3').Axes, 'title'), 'string', 'TOPAS Heterogen - TOPAS Homogen');
        
    end
end
%% Plot profiles through isoCenter

if enable(2)==1
    disp('Plotting profiles...');
    fontsize=12;
    profilex{1} = cube1(:,slicename{2},slicename{3});
    profiley{1} = permute(cube1(slicename{1},slicename{2},:),[3 2 1]);
    profilez{1} = permute(cube1(slicename{1},:,slicename{3}),[2 3 1]);
    
    profilex{2} = cube2(:,slicename{2},slicename{3});
    profiley{2} = permute(cube2(slicename{1},slicename{2},:),[3 2 1]);
    profilez{2} = permute(cube2(slicename{1},:,slicename{3}),[2 3 1]);
    
    cubes={'cube1','cube2','Dose 1','Dose 2'};
    %%
    if exist('pln','var') && ~isempty(pln)
        if strcmp(pln.bioParam.model,'none')
            l='Dose [Gy]';
        else
            l='RBE x Dose [Gy(RBE)]';
        end
    else
        l='Dose [Gy]';
    end
    
    %%
    
    hfig.profiles.fig = figure('Renderer', 'painters', 'Position', [10 50 800 800]);
    set(gcf,'Color',[1 1 1]);
    
    hfig.profiles.x = subplot(2,2,1);
    plot(profilex{1})
    hold on
    plot(profilex{2})
    xlabel('X [mm]','FontSize',fontsize)
    ylabel(l,'FontSize',fontsize);
    title('x-Profiles');
    legend({'dose w/o correction','dose with correction'},'Location','southeast')
    legend boxoff
    
    hfig.profiles.y = subplot(2,2,2);
    plot(profiley{1})
    hold on
    plot(profiley{2})
    xlabel('Y [mm]','FontSize',fontsize)
    ylabel(l,'FontSize',fontsize);
    title('y-Profiles');
    legend({'dose w/o correction','dose with correction'},'Location','southeast')
    legend boxoff
    
    hfig.profiles.z = subplot(2,2,3);
    plot(profiley{1})
    hold on
    plot(profiley{2})
    xlabel('Z [mm]','FontSize',fontsize)
    ylabel(l,'FontSize',fontsize);
    title('z-Profiles');
    legend({'dose w/o correction','dose with correction'},'Location','southeast')
    legend boxoff
    
    set(hfig.profiles.fig,'name',['Profiles:, x=',num2str(slicename{1}),'mm, y=',num2str(slicename{2}),'mm, z=',num2str(slicename{3}),'mm']);
    
end
%% Calculate and plot DVH

if enable(3)==1 && ~isempty(cst)
    disp('Calculating DVH...');
    dvh1=matRad_calcDVH(cst,cube1);
    dvh2=matRad_calcDVH(cst,cube2);
    dvh3=matRad_calcDVH(cst,cube3);
    dvhWindow = max([dvh1(1).doseGrid dvh2(1).doseGrid dvh3(1).doseGrid]);
    % Plot DVH
    disp('Plotting DVH...');
    
    cstConst = cst;
    cstMCN = cst;
    cstWED = cst;
    
    cstConst{1,2} = 'ConstRBE';
    cstMCN{1,2} = 'MCN';
    cstWED{1,2} = 'WED';
    
    hfig.dvh.fig = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
    set(gcf,'Color',[1 1 1]);
    matRad_showDVH(dvh1,cstConst,pln);
    hold on
    matRad_showDVH(dvh2,cstMCN,pln,2);
    matRad_showDVH(dvh3,cstWED,pln,3);
    xlim([0 dvhWindow*1.2])
    title('Dose Volume Histrogram, ConstRBE: solid, MCN: dotted, WED: dashed')
    
    k = 3;
    ptv = 2;
    patient = 'Boxphantom';
    mode = 'protons';
    machine = 'HIT_APM';
    doseMode = 'RBE';
    
    rT = [];
    rT{1,1} = 'patient';
    rT{1,2} = 'mean [Gy]';
    rT{1,3} = 'D_5 [Gy]';
    rT{1,4} = 'D_50 [Gy]';
    rT{1,5} = 'D_95 [Gy]';
    rT{1,6} = 'D_98 [Gy]';
    rT{1,7} = 'Delta D_95 [Gy]';
    
    qi1  = matRad_calcQualityIndicators(cstConst,pln,cube1);
    rT{k,1} = [patient,'_',mode,'_',machine,'_',doseMode,'_noCorr_',qi1(ptv).name];
    rT{k,2} = qi1(ptv).mean;
    rT{k,3} = qi1(ptv).D_5;
    rT{k,4} = qi1(ptv).D_50;
    rT{k,5} = qi1(ptv).D_95;
    rT{k,6} = qi1(ptv).D_98;
    
    qi2  = matRad_calcQualityIndicators(cstConst,pln,cube2);
    rT{k+1,1} = [patient,'_',mode,'_',machine,'_',doseMode,'_withCorr_',qi1(ptv).name];
    rT{k+1,2} = qi2(ptv).mean;
    rT{k+1,3} = qi2(ptv).D_5;
    rT{k+1,4} = qi2(ptv).D_50;
    rT{k+1,5} = qi2(ptv).D_95;
    rT{k+1,6} = qi2(ptv).D_98;
    
    qi3  = matRad_calcQualityIndicators(cstConst,pln,cube3);
    rT{k+2,1} = [patient,'_',mode,'_',machine,'_',doseMode,'_withCorr_',qi1(ptv).name];
    rT{k+2,2} = qi3(ptv).mean;
    rT{k+2,3} = qi3(ptv).D_5;
    rT{k+2,4} = qi3(ptv).D_50;
    rT{k+2,5} = qi3(ptv).D_95;
    rT{k+2,6} = qi3(ptv).D_98;
    
    rT{k,7} =   (qi2(ptv).D_95 - qi1(ptv).D_95) / qi1(ptv).D_95;
    rT{k+1,7} = (qi3(ptv).D_95 - qi1(ptv).D_95) / qi1(ptv).D_95;
    rT{k+2,7} = (qi2(ptv).D_95 - qi3(ptv).D_95) / qi3(ptv).D_95;
    
    clear qi1 qi2 qi3
    k = k + 6;
    
    disp('Done!');
    filename = 'C:\Users\homolka\Documents\#PhD\#Ergebnisse\ErgebnisseHeterogenit?tskorrektur\ModelComparison.xlsx';
    xlswrite(filename,rT,2)
end
%%

disp('Done!');

end