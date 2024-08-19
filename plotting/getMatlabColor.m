function rgb = getMatlabColor(color)

switch color
    case 'blue'
        rgb = [0, 0.4470, 0.7410];
    case 'orange'
        rgb = [0.8500, 0.3250, 0.0980];
    case 'yellow'
        rgb = [0.9290, 0.6940, 0.1250];
    case 'purple'
        rgb = [0.4940, 0.1840, 0.5560];
    case 'green'
        rgb = [0.4660, 0.6740, 0.1880];
    case 'darkGreen'
        rgb = [0.2000, 0.3098, 0.1059];
    case 'cyan'
        rgb = [0.3010, 0.7450, 0.9330];
    case 'red'
        rgb = [0.6350, 0.0780, 0.1840];
    case 'magenta'
        rgb = [1, 0, 1];
    case 'black'
        rgb = [0, 0, 0];
    case 'matradcolor'
        rgb = [0, 0.4470, 0.7410];
    case 'topascolor'
        rgb = [0.4660, 0.6740, 0.1880];
    case 'mcsquarecolor'
        rgb = [0.8500, 0.3250, 0.0980];
    case {'grey','darkGrey'}
        rgb = [0.5 0.5 0.5];
end

end

