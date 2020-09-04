% Clear the workspace and close open figure windows
close all;
clear

% Shuffle the random number generator
rng('shuffle');


%-----------------------------------
%  Possible parameters to change
%-----------------------------------

% Values used: 10, 25, 40 and 55 
trailPerStim = 10;

% Controls the spread of the individual observer data around the population
% effects
sigma = 1;

% Number of experiments
numExps = 100;

% Number of observers per condition (final sims: 5 to 15)
numObsAll = 5:15;
numObsType = length(numObsAll);

slopeEffect = linspace(0, 0.4, numObsType);

pseCD = [cd '/PSE Trials Per Stimulus Level ' num2str(trailPerStim)]; 
slopeCD = [cd '/Slope Trials Per Stimulus Level ' num2str(trailPerStim)];

if exist(pseCD, 'dir') == 0
  mkdir(pseCD)  
end

if exist(slopeCD, 'dir') == 0
  mkdir(slopeCD)  
end


%-----------------------------------
%       Setup variables
%-----------------------------------

% Sessions in the simulated experiment (X variable)
sessions = 0:5;
numSessions = length(sessions);

% Experimental condition properties
interceptEffect = 0;

% "No effect" control condition
slopeControl = 0;
interceptControl = 0;

% Generating control group line of no effect
noEffectControlLine = slopeControl .* sessions + interceptControl;





%--------------------------
%  Fitting procedure
%--------------------------

% Psychometric funtion specification
PF = @PAL_CumulativeNormal;

% Fitting parameters (fixed gamma and lambda at zero)
gamma = 0;
lambda = 0;
freeLambda = 0;

% Alpha and Beta free to vary
paramsFree = [1 1 0 0];

% Type PAL_minimize('options','help') for help
options = PAL_minimize('options');
options.TolFun = 1e-06;
options.MaxIter = 10000;
options.MaxFunEvals = 10000;

% Number of simulations
numBoots = 3;

%--------------------------
%  Fitting procedure cont
%--------------------------

% Number of steps on psychometric function (final sims: 9)
numSteps = 9;

% Range (degrees)
theRange = 10;
theHalfRange = theRange / 2;

% Setting up trials
outOfNum = ones(1, numSteps) .* trailPerStim;


%--------------------------
%  Main simulation loop
%--------------------------

c = 0;
tic

% Loop over participants
for participant = 1:length(numObsAll)
    
    % Get this number of observers in this experiment
    numObs = numObsAll(participant);
    
    % Double this for both experimental group and control
    ContT = numObs*2;
    
    % Switch point for the control group later
    ContB = ContT/2+1;
    
    % Empty matrix to filled by population sim each loop (Columns sessions,
    % rows participants (top half of rows is the experimental group, bottom
    % half of rows is the control group)
    populationPSES = nan(ContT, numSessions);
    
    %--------------------------
    %     Output matrices
    %--------------------------
    pseMat = nan(ContT, numSessions);
    slopeMat = nan(ContT, numSessions);
    
    
    
    % Loop over the magnitude of the effect
    for Mag = 1:length(slopeEffect)
        
        % Specify the experimental effect
        slope_E = slopeEffect(Mag);
        effectLine = slope_E .* sessions + interceptEffect;
        
        % Loop over the number of experiements specified (final sims: 100)
        for expNum = 1:numExps
            
            % Loop to split the the condtions for indexing to in the matrix
            % (experimental and control)
            for groupNum = 1:2
                
                % Loop over the sessions (columns of the output matrix)
                for i = 1:numSessions
                    
                    % Data for experimental group of simulated observer or
                    % the control group
                    if groupNum == 1
                        
                        % Simulating PSES with experimental effect at
                        % population (this is the effect line)
                        popEffect = effectLine;
                        
                        % Generating a variable to index off below (this is
                        % the effect for this session across all
                        % participants)
                        obsPSEs_E = randn(1, numObs) .* sigma + popEffect(i);
                        
                        % Record the effect for all participants
                        populationPSES(1:numObs, i) = obsPSEs_E;
                        
                        % Data for  control group of simulated observers
                    elseif groupNum == 2
                        
                        %simulating for PSES with no effect control
                        %at population level
                        popEffect = noEffectControlLine;
                        
                        %Generating a varaible to index off below
                        obsPSEs_C = randn(1, numObs) .* sigma + popEffect(i);
                        
                        %indexing the control condition PSES
                        populationPSES(ContB:ContT, i) = obsPSEs_C;
                        
                    end
                end
            end
            
            % For this magnitude of effect, for this experiment, we now
            % have the population data for the experiment and control
            % groups
            
            
            %------------------------------------------------
            % Now, we loop over all of the sessions
            %------------------------------------------------
            
            for sessionNum = 1:numSessions
                
                % Looping ove the Population PSes specified by column
                psesForThisSession = populationPSES(:, sessionNum);
                
                %  all Observers
                parfor thisObsNum = 1:ContT
                    
                    thisObsDone = 0;
                    
                    while thisObsDone == 0
                        
                        %--------------------------
                        %  Observer properties
                        %--------------------------
                        
                        %index of the pse for specified row of matrix above
                        thisObsPSE = psesForThisSession(thisObsNum);
                        
                        % Stimulus levels to be used in the experiment
                        stimLevels = linspace(0 - theHalfRange, 0 + theHalfRange, numSteps);
                        
                        % Setup slope variable
                        funcSigma = theRange * 0.25;
                        funcSlope = 1 / funcSigma;
                        
                        
                        %--------------------------
                        %  Simulate our observer
                        %--------------------------
                        
                        % Simulate Observer
                        simData = PAL_PF_SimulateObserverParametric([thisObsPSE, funcSlope, gamma, lambda],...
                            stimLevels, outOfNum, PF);
                        percRightward = simData ./ trailPerStim .* 100;
                        
                        
                        %--------------------------
                        %  Fit a function
                        %--------------------------
                        
                        % Parameter grid defining parameter space through which to perform a
                        % brute-force search for values to be used as initial guesses in iterative
                        % parameter search.
                        searchGrid = struct;
                        searchGrid.beta = linspace(funcSlope/4, funcSlope * 2, 20);
                        searchGrid.alpha = linspace(min(stimLevels), max(stimLevels), 20);
                        searchGrid.gamma = gamma;
                        searchGrid.lambda = lambda;
                        
                        
                        % Fit function using maximum likelihood
                        [params] = PAL_PFML_Fit(stimLevels, simData, outOfNum,...
                            searchGrid, paramsFree, PF, 'searchOptions', options);
                        
                        
                        % Check to see if this function is valid
                        [isInvalid] = validityCheck(params, stimLevels);
                        
                        if isInvalid == 0
                            thisObsDone = 1;
                        end
                        
                        
                    end

                    % Evaluating the function
                    %             extender = 2;
                    %             stimLevelsFineGrain = linspace(min(stimLevels) - extender, max(stimLevels) + extender, 100);
                    %             propCorrectModel = PF(params, stimLevelsFineGrain);
                    %             percCorrectModel = propCorrectModel .* 100;
                    %
                    
                    %------------------------------
                    %  Bootstrapping for errorbars
                    %------------------------------
                    
                    % Run the bootstrapping
                    %             [SD, paramsSim, LLSim, converged] = PAL_PFML_BootstrapParametric(...
                    %                 stimLevels, outOfNum, params, paramsFree, numBoots, PF, ...
                    %                 'searchGrid', searchGrid, 'searchOptions', options);
                    %
                    
                    % Get the PSE
                    pse = params(1);
                    pseMat(thisObsNum, sessionNum) = pse;
                    
                    % Get the PSE
                    slope = params(2);
                    slopeMat(thisObsNum, sessionNum) = slope;

                end
            end
            
            c = c+1;
            
            disp(['Participant ' num2str(participant) ' / ' num2str(length(numObsAll))...
                '| Mag ' num2str(Mag) ' / ' num2str(length(slopeEffect))...
                '| Exp ' num2str(expNum) ' / ' num2str(numExps)]);
            
            %output to specific working directory
            dlmwrite([pseCD '/PSE_ExpNum_' num2str(c) '.txt'], pseMat, '\t')
            dlmwrite([slopeCD '/Slope_ExpNum_' num2str(c) '.txt'], slopeMat, '\t')
            
        end
    end
end
toc




function [isInvalid] = validityCheck(params, stimLevels)

pse = params(1);
slope = params(2);

% See if PSEs are (1) nan, (2) inf, (3) out of range
pseLessThan = pse < min(stimLevels);
pseMoreThan = pse > max(stimLevels);
pseInf = isinf(pse);
pseNaN = isnan(pse);

% See if Slopes are (1) nan, (2) inf
slopeInf = isinf(slope);
slopeNaN = isnan(slope);

% The checking vector
checkVector = pseLessThan + pseMoreThan + pseInf + pseNaN + slopeInf + slopeNaN;
numInvalid = sum(checkVector);

% Remove the bad data
if numInvalid > 0
    isInvalid = 1;
else
    isInvalid = 0;
end

end








