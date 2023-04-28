function simulationTimeInHours = howlongdidittake(folder)

d = dir(folder);
d = d(~[d.isdir]);

startedProcessing = contains({d.name},'out','IgnoreCase',true) & ~contains({d.name},'_std','IgnoreCase',true);
startDate = [d(startedProcessing).datenum];

finishedProcessing = contains({d.name},'physicalDose.binheader','IgnoreCase',true) & ~contains({d.name},'_std','IgnoreCase',true);
finishDate = [d(finishedProcessing).datenum];

% Calculate difference and convert to hours
simulationTimeInHours = 24 * sum(finishDate - startDate);

hours = floor(simulationTimeInHours);
minutes = round((simulationTimeInHours - hours)*60);

disp(['Simulation time was ' num2str(hours), 'h:' num2str(minutes) 'm with ' num2str(sum(finishedProcessing)) ' batches.'])

end

