%% distributeDataForSurrogate
%
% Auxiliary routine that receives input data from DOE and process it
% accordingly to surrogate building
%
%% Inputs
% The input is a structure containing the following fields:
%%
%
% * doeBuild - this the design of experiment data (design + observed data)
% * designIndex - column indexes of design sites in doeBuild
% * fobjIndex - column index of objective function data in doeBuild
% * numConst - number of constraints in the model
% * constIndex - column indexes of constraint data in doeBuild
% * constLimit - Limit value for each constraint like $g(x) \leq
% constLimit$
% * constType - each constraint has to be classified in one the following
% types: 
% 
% # constType = 0 for constraints of $g(x) \leq constLimit0$
% # constType = 1 for constraints of $g(x) \geq constLimit$
% # constType = 2 for constraints of $lowerConstLimit \leq g(x) \leq  upperConstLimit$
%
%% Outputs
% The output arguments are:
%
% * xlhs - design experiment data. Same as doeBuild(:,designIndex)
% * fobs - objective function experiment data. Same as
% doeBuild(:,fobjIndex)
% * gobs - constraint experiment data.

function [xlhs,fobs,gobs] = distributeDataForSurrogate(dataForSurrogate)
 
doeBuild = dataForSurrogate.doeBuild;
designIndex = dataForSurrogate.designIndex;
fobjIndex = dataForSurrogate.fobjIndex;
numConst = dataForSurrogate.numConst;
constIndex = dataForSurrogate.constIndex;
constType = dataForSurrogate.constType;
constLimit = dataForSurrogate.constLimit;

% Seperate data into proper variables for surrogate building
xlhs = doeBuild(:,designIndex);
fobs = doeBuild(:,fobjIndex);
constObs = doeBuild(:,constIndex);

% get index of constraint type
type2ConstIndex = find(constType == 2);

% check if there are type 2 constraints and do some preparation
if ~isempty(type2ConstIndex) % type 2 constraint specified
    % Duplicate the 2s values
    constType = insertrows(constType',constType(type2ConstIndex)',type2ConstIndex)';
    constObs  = insertrows(constObs',constObs(:,type2ConstIndex)',type2ConstIndex)';
    if numel(constType) ~= numel(constLimit)
        error('%d constraints were specified, %d of them being type 2 (interval) constraint. However, %d value limits for constraints were specified. For each type 2 constraint, an additional constraint limit value is needed. In this case, %d limits are need. Check the number of limits or types of constraints specified.',numConst,numel(type2ConstIndex),length(constLimit),numConst+numel(type2ConstIndex));
    end
else
    if numel(constType) ~= numel(constLimit) || numConst ~= numel(constType) || numConst ~= numel(constLimit)
        error('%d constraints were specified. However, %d value limits for constraints and %d types of constraints were specified. They must be the same.',numConst,numel(constLimit),numel(constType));
    end
end

constData = NaN(size(constObs)); % Initialization
type0ConstIndex = find(constType == 0);
type1ConstIndex = find(constType == 1);
type2ConstIndex = find(rem(find(constType == 2),2)); % get the first position of 2 in the pairs

% Distribute the constraint data
if ~isempty(type0ConstIndex)
    constData(:,type0ConstIndex) = constObs(:,type0ConstIndex) - constLimit(type0ConstIndex);
end

if ~isempty(type1ConstIndex)
    constData(:,type1ConstIndex) = constLimit(type1ConstIndex) - constObs(:,type1ConstIndex);
end

if ~isempty(type2ConstIndex)
    constData(:,type2ConstIndex) = constObs(:,type2ConstIndex) - constLimit(type2ConstIndex + 1);
    constData(:,type2ConstIndex + 1) = constLimit(type2ConstIndex) - constObs(:,type2ConstIndex);
end

gobs = constData;   % set constraints to c(x) <= 0
end