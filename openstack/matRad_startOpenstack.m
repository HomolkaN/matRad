function matRad_startOpenstack(folderToStart)
[usr,pw,ip] = getOpenStackCredentials();

call = strjoin(['wsl plink -pw ''',pw,''' ',usr,'@',ip,' ''cd files; ./calcAll.sh ', folderToStart,''''],'');
[~,errMsg] = system(call);

if ~isempty(strfind(lower(errMsg),'not exist'))
    error(convertCharsToStrings(errMsg(1:end-11)));
elseif ~isempty(strfind(lower(errMsg),lower('All done!')))
    disp('Successfully started TOPAS.'); 
end
end

