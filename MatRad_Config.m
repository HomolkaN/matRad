classdef MatRad_Config < handle
    % MatRad_Config MatRad Configuration class
    % This class is used globally through Matlab to handle default values and
    % logging and is declared as global matRad_cfg.
    % Usage:
    %    matRad_cfg = MatRad_Config.instance();
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        
        %Logging
        logLevel = 3; %1 = only Errors, 2 = with Warnings, 3 = Info output, 4 = deprecation warnings, 5 = debug information
        keepLog = false; %Stores the full log in memory
        writeLog = false; %Writes the log to a file on-the-fly

        %Default Properties
        propDoseCalc;
        propOpt;
        propMC;
        propStf;
        propHeterogeneity;

        defaults;

        %Disable GUI
        disableGUI = false;
    end

    properties (SetAccess = private)
        messageLog = {};
        logFileHandle;

        %For storing the Environment & its version
        env;
        envVersion;
        isOctave; %Helper bool to check for Octave
        isMatlab; %Helper bool to check for Matlab
        matRad_version; %MatRad version string
    end

    properties (SetAccess = private)
        matRadRoot;
    end
    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Constructs an instance of this class.
            %  The configuration is implemented as a singleton and used globally
            %  Therefore its constructor is private
            %  For instantiation, use the static MatRad_Config.instance();

            obj.matRadRoot = fileparts(mfilename('fullpath'));
            addpath(genpath(obj.matRadRoot));

            %Set Version
            obj.getEnvironment();
            obj.matRad_version = matRad_version();
            
            %Just to catch people messing with the properties in the file
            if ~isempty(obj.writeLog) && obj.writeLog
                logFile = [obj.matRadRoot filesep 'matRad.log'];
                obj.logFileHandle = fopen(logFile,'a');
            end

            %Call the reset function for remaining inatialization
            obj.reset();
        end

        function displayToConsole(obj,type,formatSpec,varargin)
            %displayToConsole lowest-level logging function for matRad.
            %   Display to console will be called from the public wrapper
            %   functions dispError, dispWarning, dispInfo, dispDebug
            %
            %  input
            %    type:			type of the log information.
            %                   Needs to be one of 'error', 'warning', 'info' or 'debug'.
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec

            if nargin < 4
                forwardArgs = {formatSpec};
            else
                forwardArgs = [{formatSpec},varargin(:)'];
            end

            if obj.keepLog
                obj.messageLog{end+1,1} = upper(type);
                obj.messageLog{end,2} = sprintf(forwardArgs{:});
            end

            switch type
                case{'info'}
                    if obj.logLevel >= 3
                        fprintf(forwardArgs{:});
                    end
                case{'debug'}
                    if obj.logLevel >= 5
                        forwardArgs{1} = ['DEBUG: ' forwardArgs{1}];
                        fprintf(forwardArgs{:});
                    end
                case{'dep'}
                    if obj.logLevel >= 4
                        forwardArgs{1} = ['DEPRECATION WARNING: ' forwardArgs{1}];
                        warning(forwardArgs{:});
                    end
                case{'warning'}
                    if obj.logLevel >= 2
                        warning(forwardArgs{:});
                    end
                case {'error'}
                    if obj.logLevel >= 1
                        %We create an error structure to later clean the
                        %stack trace from the last two files/lines (i.e.,
                        %this function / file)

                        err.message = sprintf(forwardArgs{:});
                        err.identifier = 'matRad:Error';
                        err.stack = dbstack(2);
                        error(err);

                    end
                otherwise
                    error('Log type %s not defined!',type);
            end

            if obj.writeLog
                fprintf(obj.logFileHandle,forwardArgs{:});
            end
        end 
    end

    methods
        function reset(obj)
            %Set all default properties for matRad's computations
            obj.setDefaultProperties();
        end

        function setDefaultProperties(obj)
            %setDefaultProperties set matRad's default computation
            %   properties
            %  input

            obj.propStf.defaultLongitudinalSpotSpacing = 2;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding

            obj.propDoseCalc.defaultResolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.propDoseCalc.defaultLateralCutOff = 0.995; %[rel.]
            obj.propDoseCalc.defaultGeometricCutOff = 50; %[mm]
            obj.propDoseCalc.defaultKernelCutOff = Inf; %[mm]
            obj.propDoseCalc.defaultSsdDensityThreshold = 0.05; %[rel.]
            obj.propDoseCalc.defaultUseGivenEqDensityCube = false; %Use the given density cube ct.cube and omit conversion from cubeHU.
            obj.propDoseCalc.defaultIgnoreOutsideDensities = true; %Ignore densities outside of cst contours
            obj.propDoseCalc.defaultVoxelSubIx = []; %Allows specification of a subindex list for dose calculation, empty by default means automatic setting
            obj.propDoseCalc.defaultUseCustomPrimaryPhotonFluence = false; %Use a custom primary photon fluence
            obj.propDoseCalc.defaultCalcLET = true; %calculate LETs for particles

            obj.propDoseCalc.defaultAirOffsetCorrection = true;

            obj.propOpt.defaultMaxIter = 500;
            obj.propOpt.defaultRunDAO = 0;
            obj.propOpt.defaultRunSequencing = 0;

            obj.propMC.defaultCarbonEnergySpread = 0; %[%]

            obj.propMC.ompMC_defaultHistories = 1e6;
            obj.propMC.ompMC_defaultOutputVariance = false;

            obj.propMC.direct_defaultHistories = 1e7;
            obj.propMC.particles_defaultHistories = 2e4;
            obj.propMC.MCsquare_defaultHistories = 1e8;

            %obj.propMC.default_photon_engine = 'ompMC';
            obj.propMC.default_proton_engine = 'MCsquare';
            obj.propMC.default_carbon_engine = 'TOPAS';

            % Default settings for TOPAS
            obj.propMC.defaultExternalCalculation = false;
            obj.propMC.topas_defaultNumBatches = 5;
            obj.propMC.defaultCalcDij = false;

            obj.propHeterogeneity.defaultType = 'complete'; % 'depthBased','voxelwise'
            obj.propHeterogeneity.defaultCalcHetero = true;
            obj.propHeterogeneity.defaultUseOriginalDepths = false; % use original base data depths instead of fitted ones
            obj.propHeterogeneity.defaultModulateBioDose = true; % directly modulate alpha beta curves for RBE
            obj.propHeterogeneity.defaultModulateLET = false;

            obj.propHeterogeneity.defaultModPower = 800;
            obj.propHeterogeneity.sampling.defaultMode = 'matRad';
            obj.propHeterogeneity.sampling.defaultMethod = 'binomial';
            obj.propHeterogeneity.sampling.defaultUseContinuous = true;
            obj.propHeterogeneity.sampling.defaultNumOfSamples = 50;
            obj.propHeterogeneity.sampling.defaultContinuous = true;

            
            % default properties for fine sampling calculation
            obj.propDoseCalc.fineSamplingProperties.defaultSigmaSub = 1;
            obj.propDoseCalc.fineSamplingProperties.defaultN = 21;
            obj.propDoseCalc.fineSamplingProperties.defaultMethod = 'russo';

            obj.disableGUI = false;
            
            obj.defaults.samplingScenarios = 25;
        end

        %%For testing
        function setDefaultPropertiesForTesting(obj)
            %setDefaultPropertiesForTesting sets matRad's default
            %properties during testing to reduce computational load

            obj.logLevel   = 3; %Omit output except errors

            obj.propStf.defaultLongitudinalSpotSpacing = 20;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding

            obj.propDoseCalc.defaultResolution = struct('x',5,'y',6,'z',7); %[mm]
            obj.propDoseCalc.defaultGeometricCutOff = 20;
            obj.propDoseCalc.defaultLateralCutOff = 0.8;
            obj.propDoseCalc.defaultKernelCutOff = 20; %[mm]
            obj.propDoseCalc.defaultSsdDensityThreshold = 0.05;
            obj.propDoseCalc.defaultUseGivenEqDensityCube = false; %Use the given density cube ct.cube and omit conversion from cubeHU.
            obj.propDoseCalc.defaultIgnoreOutsideDensities = true;
            obj.propDoseCalc.defaultVoxelSubIx = []; %Allows specification of a subindex list for dose calculation, empty by default means automatic setting
            obj.propDoseCalc.defaultUseCustomPrimaryPhotonFluence = false; %Use a custom primary photon fluence
            obj.propDoseCalc.defaultCalcLET = true; %calculate LET for particles
            
            % default properties for fine sampling calculation
            obj.propDoseCalc.fineSamplingProperties.sigmaSub = 2;
            obj.propDoseCalc.fineSamplingProperties.N = 5;
            obj.propDoseCalc.fineSamplingProperties.method = 'russo';
            
            obj.propOpt.defaultMaxIter = 10;

            obj.propMC.ompMC_defaultHistories = 100;
            obj.propMC.ompMC_defaultOutputVariance = true;

            obj.propMC.particles_defaultHistories = 100;
            obj.propMC.direct_defaultHistories = 100;

            %obj.propMC.default_photon_engine = 'ompMC';
            obj.propMC.default_proton_engine = 'MCsquare';
            obj.propMC.default_carbon_engine = 'TOPAS';

            % Default settings for TOPAS
            obj.propMC.topas_defaultNumBatches = 5;

            obj.propHeterogeneity.defaultType = 'complete'; % 'depthBased','voxelwise'
            obj.propHeterogeneity.defaultCalcHetero = true;
            obj.propHeterogeneity.defaultUseOriginalDepths = true; % use original base data depths instead of fitted ones
            obj.propHeterogeneity.defaultUodulateBioDose = true; % use alpha beta curves for RBE

            obj.defaults.samplingScenarios = 2;

            obj.disableGUI = true;
        end

        function dispDebug(obj,formatSpec,varargin)
            %dispDebug print debug messages (log level >= 4)
            %  input
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec

            obj.displayToConsole('debug',formatSpec,varargin{:});
        end

        function dispInfo(obj,formatSpec,varargin)
            %dispInfo print information console output (log level >= 3)
            %  input
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec
            obj.displayToConsole('info',formatSpec,varargin{:});
        end

        function dispError(obj,formatSpec,varargin)
            %dispError print errors (forwarded to "error" that will stop the program) (log level >= 1)
            %  input
            %    formatSpec: 	string to print using format specifications
            %                   similar to 'error'
            %    varargin:   	variables according to formatSpec
            obj.displayToConsole('error',formatSpec,varargin{:});
        end

        function dispWarning(obj,formatSpec,varargin)
            %dispError print warning (forwarded to 'warning') (log level >= 2)
            %  input
            %    formatSpec: 	string to print using format specifications
            %                   similar to 'warning'
            %    varargin:   	variables according to formatSpec
            obj.displayToConsole('warning',formatSpec,varargin{:});
        end

        function dispDeprecationWarning(obj,formatSpec,varargin)
            %dispDeprecationWarning wrapper for deprecation warnings forwarded to displayToConsole
            obj.displayToConsole('dep',formatSpec,varargin{:});
        end

        function obj = writeLogToFile(obj,filename)
            %writeLogToFile writes the log kept in MatRad_Config to file.
            %  Note that the switch keepLog must be enabled for MatRad_Config to store all logging output.

            singleString = '%s: %s\n';
            fID = fopen(filename,'w');
            fprintf(fID,repmat(singleString,1,size(obj.messageLog,1)),obj.messageLog{:});
            fclose(fID);
        end

        function set.logLevel(obj,newLogLevel)
            %%Property set methods for logLevel
            minLevel = 1;
            maxLevel = 5;
            if newLogLevel >= minLevel && newLogLevel <= maxLevel
                obj.logLevel = newLogLevel;
            else
                obj.dispError('Invalid log level. Value must be between %d and %d',minLevel,maxLevel);
            end
        end

        function set.writeLog(obj,writeLog)
            if writeLog
                logFile = [obj.matRadRoot filesep 'matRad.log'];
                obj.logFileHandle = fopen(logFile,'a');
                obj.writeLog = true;
            else
                fclose(obj.logFileHandle);
                obj.writeLog = false;
            end
        end

        function getEnvironment(obj)
            % getEnvironment function to get the software environment
            %   matRad is running on

            obj.isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
            obj.isMatlab = ~obj.isOctave;

            if obj.isOctave
                obj.env = 'OCTAVE';
                obj.envVersion = OCTAVE_VERSION;
            else
                obj.env = 'MATLAB';
                vData = ver(obj.env);
                obj.envVersion = vData.Version;

            end
        end

        function pln = loadDefaultParam(obj,pln)
            % function to load all non-set parameters into pln struct
            structs = {'propHeterogeneity','propDoseCalc','propOpt'};
            for propName = structs
                if ~isfield(pln,propName{1})
                    pln.(propName{1}) = struct();
                end

                fnames = fieldnames(obj.(propName{1}));
                for f = 1:length(fnames)
                    if contains(fnames{f},'default')
                        cutName = [lower(fnames{f}(8)) fnames{f}(9:end)];
                        if ~isfield(pln.(propName{1}),cutName)
                            pln.(propName{1}).(cutName) = obj.(propName{1}).(fnames{f});
                        end
                    else
                        if ~isfield(pln.(propName{1}),fnames{f})
                            pln.(propName{1}).(fnames{f}) = struct();
                        end
                        subfields = fieldnames(obj.(propName{1}).(fnames{f}));
                        for s = 1:length(subfields)
                            if contains(subfields{s},'default')
                                if length(subfields{s})==8
                                    cutName = [subfields{s}(8)];
                                else
                                    cutName = [lower(subfields{s}(8)) subfields{s}(9:end)];
                                end
                                if ~isfield(pln.(propName{1}).(fnames{f}),cutName)
                                    pln.(propName{1}).(fnames{f}).(cutName) = obj.(propName{1}).(fnames{f}).(subfields{s});
                                end
                            end
                        end
                    end
                end

            end

            % Turn off modulateBioDose in case no suitable RBE model is used
            if isfield(pln,'bioParam') && ~pln.bioParam.bioOpt
                pln.propHeterogeneity.modulateBioDose = 0;
            end

            % Monte Carlo settings
            if strcmp(pln.propHeterogeneity.sampling.mode,'TOPAS')
                pln.propMC.engine = 'TOPAS';
            end
            if isfield(pln,'propMC')
                if isfield(pln.propMC,'outputVariance')
                    obj.dispWarning('Variance scoring for TOPAS not yet supported.');
                end

                if strcmp(pln.radiationMode,'protons')
                    engines = {'TOPAS','MCsquare'};
                    if ~isfield(pln.propMC,'engine') || ~any(strcmp(pln.propMC.engine,engines))
                        obj.dispInfo('Using default proton MC engine "%s"\n',obj.propMC.default_proton_engine);
                        pln.propMC.engine = obj.propMC.default_proton_engine;
                    end


                elseif strcmp(pln.radiationMode,'carbon') || strcmp(pln.radiationMode,'helium')
                    if ~isfield(pln.propMC,'engine')
                        pln.propMC.engine = obj.propMC.default_carbon_engine;
                    end

                else
                    obj.dispError('MC only implemented for protons, helium and carbon ions (only TOPAS).');
                end
            end
            switch pln.propMC.engine
                % number of histories per beamlet (nCasePerBixel > 1),
                % max stat uncertainity (0 < nCasePerBixel < 1)
                % set number of particles simulated per pencil beam
                case 'MCsquare'
                    if ~isfield(pln.propMC,'histories') || (pln.propMC.histories==obj.propMC.MCsquare_defaultHistories)
                        pln.propMC.histories = obj.propMC.MCsquare_defaultHistories;
                        obj.dispInfo('Using default number of Histories per Bixel: %d\n',pln.propMC.histories);
                    end
                case 'TOPAS'
                    if ~isfield(pln.propMC,'externalCalculation')
                        pln.propMC.externalCalculation = obj.propMC.defaultExternalCalculation;
                    end
                    if ~isfield(pln.propMC,'histories') || (pln.propMC.histories==obj.propMC.particles_defaultHistories)
                        pln.propMC.histories = obj.propMC.particles_defaultHistories;
                        obj.dispInfo('Using default number of Histories per Bixel: %d\n',pln.propMC.histories);
                    end
            end
        end
    end

    methods(Static)

        function obj = instance()
            %instance creates a singleton instance of MatRad_Config
            %  In MatRad_Config, the constructor is private to make sure only on global instance exists.
            %  Call this static functino to get or create an instance of the matRad configuration class
            persistent uniqueInstance;

            if isempty(uniqueInstance)
                obj = MatRad_Config();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
               
        function obj = loadobj(sobj)
        % Overload the loadobj function to allow downward compatibility
        % with workspaces which where saved as an older version of this class
        
            function basic_struct = mergeStructs(basic_struct, changed_struct)
                % nested function for merging the properties of the loaded
                % obj into a new obj.
                % Merges two structs, including nestes structs, by overwriting 
                % the properties of basic_struct with the changed properties in changed_struct
                fields = fieldnames(basic_struct);
                for k = 1:length(fields)  
                    disp(fields{k});
                    if(isfield(changed_struct, fields{k}))                 
                        if isstruct(changed_struct.(fields{k})) && isstruct(basic_struct.(fields{i}))        
                            basic_struct.(fields{k}) = mergeStructs(basic_struct.(fields{k}), changed_struct.(fields{i}));
                        else
                            basic_struct.(fields{k}) = changed_struct.(fields{k});
                        end
                    end
                end
            end
            
            % If the saved object is loaded as a struct there was a problem
            % with the generic loading process most likly a version-conflict
            % regarding the structs, in order to fix this, do a custom 
            % loading process including recursivly copying the conflicting structs 
            if isstruct(sobj)
                warning('The  loaded object differs from the current MatRad_Config class, resuming the loading process with the overloaded loadobj function!');
                obj = MatRad_Config(); 
                % Use a metaclass object to get the properties because
                % Octave <= 5.2 doesn't have a properties function
                props = {metaclass(obj).PropertyList.Name};
                % Throw warning if the version differs and remove the
                % matRad_version field from the loaded struct, in order to
                % not overwrite the version later
                if (isfield(sobj, 'matRad_version') && ~(strcmp(obj.matRad_version, sobj.matRad_version)))
                    warning('MatRad version or git Branch of the loaded object differs from the curret version!');
                    sobj = rmfield(sobj, 'matRad_version');
                end
                % Itterate over the properties of the newly created MatRad_Config object
                for i = 1:length(props)
                    % check if the field exists in the loaded object
                    if(isfield(sobj,props{i}))
                        objField = obj.(props{i});
                        sobjField = sobj.(props{i});
                        % If field from loaded object and from the newly
                        % created object are equal skip it, else copy the
                        % value of the loaded object and if it's a struct
                        % check it's field recursively
                        if ~(isequal(sobjField, objField))
                            if (isstruct(sobjField) && isstruct(objField))
                                retStruct = mergeStructs(objField,sobjField);
                                obj.(props{i}) = retStruct;
                            else
                                obj.(props{i}) = sobjField;
                            end
                        end
                    end
                end
            else
                obj = sobj;
            end     
        end
        
        
    end
end

