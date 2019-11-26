function sendMCtoOpenStack(project_name)% This script saves a workspace (with optimized dose distribution!), copies
% it to an OpenStack TOPAS MC machine and starts the simulation
% Requirements:
% - Matlab SSH client: https://de.mathworks.com/matlabcentral/fileexchange/35409-ssh-sftp-scp-for-matlab-v2
% - At best access to the E071-FSE-McMegaMuscle Project in the DKFZ OpenStack
% - A running OpenStack machine (I suggest the flavor dkfz-28.64 --> 28 cores, 64GB RAM) built from the image CentOS_DKFZ_MC_TOPAS (image ID = 227dd0f2-9acb-412f-bd82-f3b6deb29540)
% - An ssh key pair to connect to this machine (standard linux, not putty format)
% - The floating IP for this machine (= hostname)


folder_name = ['matfiles/Boxphantom_protons/' project_name]; %Project name (matches the folder name that is created on the host machine)
hostname = '10.128.129.131'; %OpenStack host in project E071-FSE-McMegaMuscle which uses the image CentOS_DKFZ_MC_TOPAS
user = 'centos';
key_private = 'D:\homolka\Documents\#PhD\openstack\keys\id_rsa'; %Local path to private key file (standard linux rsa, NOT in putty format)
key_public = 'D:\homolka\Documents\#PhD\openstack\keys\id_rsa.pub'; %Local path to public key file

fracHistories = 1e-7; %Fraction of particle histories to simulate
numThreads = 28; %Number of threads to be used (28 is optimal for the dkfz-28.64 flavor used

matFileName = [folder_name '.mat'];

%% Save workspace
save(matFileName,'-v7'); % Needs octave compatible format

%% Establish connection
ssh2_conn = ssh2_config_publickey(hostname,user,key_private,key_public);

%% Create Folder
%Check for existing folder
[ssh2_conn,command] = ssh2_command(ssh2_conn,'ls');
if ismember(folder_name,command)
    error('Project already exists on the MC machine! Choose different Name!');
end
ssh2_conn = ssh2_command(ssh2_conn, ['mkdir ' folder_name]);

%% File Copy
[ssh2_conn] = scp_put(ssh2_conn,matFileName,['~/' folder_name '/']);


%% Execute
%Forwards the MC simulation and waits for the results, console output is
%stored in "output"
execCommand = ['cd ' folder_name ' && source ~/topas/setup_env.sh && forwardMC ' matFileName ' --fracHistories ' num2str(fracHistories) ' --threads ' num2str(numThreads) ' --label ' folder_name];
[ssh2_conn,output] = ssh2_command(ssh2_conn, execCommand);

%ssh2_conn = ssh2_command(ssh2_conn,['nohup forwardMC ' matFileName ' --fracHistories ' num2str(fracHistories) ' --threads ' num2str(numThreads) ' &']);

%% Retreive Result

%TOPAS saves into weird folders using a linux PID. We need to find this
%folder from the console output via a regular expression on the log
[tokens,match] = regexp(output{5},'PID ([0-9]+)','tokens','match');
PID = tokens{1}{1};

%Obtain the mat file with the monte carlo results
[ssh2_conn] = scp_get(ssh2_conn,'workspace_withMCresults.mat',[],['~/' folder_name '/data_PID' PID '_' folder_name '/']);

%% Close the connection
ssh2_conn = ssh2_close(ssh2_conn);

%% Load Data
clear;

load 'workspace_withMCresults.mat';

matRadGUI;
end