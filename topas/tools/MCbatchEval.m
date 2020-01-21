folder = dir('files/data_P*');
MCresults = cell(length(folder),1);
MCresults{length(folder),1} = [];
for i = 1:length(folder)
    tmp = load(['files/' folder(i).name '/MCdata.mat'],'resultMC');
    MCresults{i} = tmp.resultMC.MC_physicalDose;
    clear tmp
end
%% Average all MC dose cubes
arrayDim = size(MCresults{1});
sumDoseMC = zeros(arrayDim(1),arrayDim(2),arrayDim(3)) ;
for i = 1:length(folder)
    sumDoseMC = sumDoseMC + MCresults{i};
end
meanPhysDoseMC = sumDoseMC / length(folder);
% MC.heterogen.physicalDose = meanPhysDoseMC;
TOPAS = meanPhysDoseMC;
