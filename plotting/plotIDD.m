function plotIDD(varargin)
%PLOTIDD Summary of this function goes here
%   Detailed explanation goes here
% figure
hold on

profile = 0;
linestyles1 = {'-',':','-.'};
linestyles2 = {getMatlabColor('blue'),getMatlabColor('orange'),getMatlabColor('green')};


if sum(cellfun(@(x) iscell(x), varargin)) == 2
    counter = length(varargin)-2;
elseif sum(cellfun(@(x) iscell(x), varargin)) == 0
    counter = length(varargin);
    varargin{end+1} = repmat(linestyles1,1,round(counter/2));
    varargin{end+1} = repmat(linestyles2,1,round(counter/2));
end
% for j = 1:length(varargin)
%     if isstruct(varargin{j})
%         counter = counter + 1;
%     end
% end

for i = 1:counter
    plot(matRad_calcIDD(varargin{i},profile),'LineWidth',1.5, 'DisplayName', inputname(i),'LineStyle',varargin{counter+1}{i},'Color',varargin{counter+2}{i})
end
if profile
    % ylabel('profile / RBExD (Gy)')
    ylabel('profile / physical dose (Gy)')
else
    % ylabel('IDD / RBExD (Gy)')
    ylabel('IDD / physical dose (Gy)')
end
xlabel('depth [voxel]')

% xlim([0 size(varargin{i}.physicalDose,1)*0.8])
legend('Interpreter', 'none','Location','northwest')
end

