function sortBasedata(pathToBasedata)

load(pathToBasedata,'machine')

fnames = lower(fieldnames(machine.data));

alphaBeta = any(cellfun(@(teststr) ~isempty(strfind(fnames,teststr)), {'alpha','beta','LET','AbsNb'}));
fnamesAB = fnames;
fnamesAB(alphaBeta) = deal({'zzzz'});
[~,idx] = sort(fnamesAB);

machine.data = orderfields(machine.data, idx);

save(pathToBasedata,'machine')

end