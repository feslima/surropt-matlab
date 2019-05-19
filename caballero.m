function [xlhs, gobs, fobs, xstar, fstar] = caballero(sampleFunction, prob_struct)
% TODO: prob_struct field verification (sanitation)
doeBuild = prob_struct.doeBuild;
regModel = prob_struct.regModel;
corModel = prob_struct.corModel;

[nInitialPoints, ~] = size(doeBuild);

lb = prob_struct.lb;
ub = prob_struct.ub;

designIndex = prob_struct.designIndex;
fobjIndex = prob_struct.fobjIndex;
constIndex = prob_struct.constIndex;

numConst = numel(constIndex);
constLimit = prob_struct.constLimit;

conTol = prob_struct.conTol;
constType = prob_struct.constType;

tol1 = prob_struct.tol1;
tol2 = prob_struct.tol2;
tolContraction = prob_struct.tolContraction;

firstContractionFactor = prob_struct.firstContractionFactor;
secondContractionFactor = prob_struct.secondContractionFactor;
MaxFunEvals = prob_struct.MaxFunEvals;
penaltyFactor = prob_struct.penaltyFactor;

nlp_solver = prob_struct.nlp_solver;

% load auxiliary functions and modules
addpath('..\dace')
addpath('..\lhs')
addpath('..\sqp')
addpath('..\ipopt')

% Encase the data needed for surrogate building
dataForSurrogate.doeBuild = doeBuild;
dataForSurrogate.designIndex = designIndex;
dataForSurrogate.fobjIndex = fobjIndex;
dataForSurrogate.numConst = numConst;
dataForSurrogate.constIndex = constIndex;
dataForSurrogate.constType = constType;
dataForSurrogate.constLimit = constLimit;

% shape the data
[xlhs,fobs,gobs] = distributeDataForSurrogate(dataForSurrogate);

% Construct the surrogate - calibrate kriging parameters
[krmodelfobj,krmodelcon] = buildSurrogate(xlhs,[fobs gobs],regModel,corModel);

% iteration count
k = 1; j = 1;

% Search for feasible points in the initial sampling
bestSampledIndex(k,1) = find(min(fobs(all(gobs <= conTol, 2))) == fobs); % upkeep
if isempty(bestSampledIndex)
    error('Search for a initial feasible point. Implement ISC 14 from Sasena.');
end

x0 = xlhs(bestSampledIndex,:);

ubopt = ub; lbopt = lb; hlb = lb; hub = ub; dlb = lb; dub = ub;
moveNum = 0; contractNum = 0;
funEvals = 0;
fprintf('k\tj\t xjk \n');
tic

largeContractPerformed = false;

while true
    [xjk, fjk, exitflag] = optimizeNLP(krmodelfobj, krmodelcon, x0, ...
        lbopt, ubopt, nlp_solver);
    
    if exitflag < 0
        fprintf(2,'%d\t%d\t | %s| (infeasible)\n',k, j, sprintf('%8.4f ',xjk))
    else
        fprintf('%d\t%d\t | %s| (feasible)\n',k, j, sprintf('%8.4f ',xjk))
    end
    
    % Check for maximum number of function evaluations
    if funEvals >= MaxFunEvals
        warning('Maximum number of function evaluations achieved!');
        if ~isempty(fobs(all(gobs <= conTol, 2)))
            feasibleIndex = find(min(fobs(all(gobs <= conTol, 2))) == fobs);
            fprintf('\nBest feasible value found: %8.6f\nat point\n',fobs(feasibleIndex))
            disp(xlhs(feasibleIndex,:));
            idx = rangesearch(xlhs,xlhs(feasibleIndex,:),0.01*norm(lb-ub));
            fprintf('%d of points are within 1%% range of this point\n',numel(idx{:})-1);
            fprintf('Number of funcion evaluations: %d\n',funEvals);
        end
        toc
        break
    end
    
    if ~ismember(xjk,xlhs,'rows') % if the solution found is not repeated
        % sample in the point (xj)k and update kriging (no-reoptimize)
        [xSampled,fobsSampled,gobsSampled] = sampleModel(xjk,sampleFunction, dataForSurrogate);
        xlhs(end+1,:) = xSampled;
        fobs(end+1) = fobsSampled;
        gobs(end+1,:) = gobsSampled;
        
    else
        fprintf(2', 'Optimizer couldn''t improve from last iteration.')
        break
    end
    
    funEvals = funEvals + 1;
    
    if k~=1 && j == 1
        [krmodelfobj,krmodelcon] = buildSurrogate(xlhs,[fobs gobs],...
            regModel,corModel);
    else % do not update kriging parameters, just insert points
        [krmodelfobj,krmodelcon] = buildSurrogate(xlhs,[fobs gobs],...
            regModel,corModel,krmodelfobj,krmodelcon);
    end
    
    % Starting from xjk solve the NLP to get xj1k
    [xj1k, fj1k, exitflag] = optimizeNLP(krmodelfobj, krmodelcon, xjk, ...
        lbopt, ubopt, nlp_solver);
    
    if norm(xjk - xj1k)/norm(lb-ub) >= tol1 % did improve, keep going
        xjk = xj1k;
        j = j + 1;
    else
        if k == 1
            % initialize xstar with best incumbent solution
            xstar(k, :) = x0;
            gstar(k, :) = gobs(bestSampledIndex,:);
            fstar(k, 1) = fobs(bestSampledIndex);
            
            moveTrack(k,1) = -1;  % -1 for first iteration
        end
        
        % find best value inserted
        xlhsInserted = xlhs((nInitialPoints+1):end,:);
        fobsInserted = fobs((nInitialPoints+1):end,:);
        gobsInserted = gobs((nInitialPoints+1):end,:);
        
        bestInserted = find(min(fobsInserted(all(gobsInserted <= conTol, 2))) == fobsInserted);
        
        if ~isempty(bestInserted)
            if numel(bestInserted) > 1
                bestInserted = bestInserted(end);
            end
            
            % TODO # 1:
            % "Eventually if in two consecutives iterations we get the same
            % values of the independent variables with no improvement, we
            % cannot guarantee an optimum in the actual model if the
            % predicted errors in the surrogate model are greater than the
            % tolerance."
            %
            % "In noise systems it is not possible to verify if the
            % gradient of the metamodel matches the gradient of the true
            % function. In this case, the stopping criterion is based on
            % the assumption that if two sucessive major iterations (AT
            % LEAST ONE CONTRACTION MUST BE PERFORMED) the optimal solution
            % is the same, we expect that the gradients also match the
            % 'true gradients'."
            
            xstar(k + 1,:) = xlhsInserted(bestInserted,:);
            gstar(k + 1,:) = gobsInserted(bestInserted,:);
            fstar(k + 1,1) = fobsInserted(bestInserted,:);
        else
            % if there is no best sampled, just insert the last value found
            xstar(k + 1,:) = xlhs(end,:);
            gstar(k + 1,:) = gobs(end,:);
            fstar(k + 1,1) = fobs(end,:);
        end
        
        % Select the best point to be centered
        bestFeasibleStar = find(min(fstar(all(gstar <= conTol, 2))) == fstar);
        
        if numel(bestFeasibleStar) > 1 % if the find returns multiple values, get the last
            bestFeasibleStar = bestFeasibleStar(end);
        end
        
        xstark = xstar(bestFeasibleStar,:);
        
        % check for hypercube limits
        if contractNum == 0
            contractFactor = firstContractionFactor;
        else
            contractFactor = secondContractionFactor;
        end
        
        [hub, hlb, exitflagref,ubopt,lbopt] = refineHypercube(xstark,dlb,hlb,dub,hub,contractFactor,tolContraction);
        
        if exitflagref == 2 || exitflagref == 4
            fprintf(2,'--------------- Move! ---------------\n')
            moveNum = moveNum + 1;
            moveTrack(k + 1,1) = 0; % 0 for hypercube movement
        end
        
        if exitflagref == 1 || exitflagref == 3
            fprintf(2,'------------- Contract! -------------\n')
            contractNum = contractNum + 1;
            moveTrack(k + 1,1) = 1; % 1 for hypercube contraction
            % resample the reduced hypercube
            [samples,fobsSampled,gobsSampled] = getSamples(xlhsInserted,fobsInserted,gobsInserted,hlb,hub);
            
            % each contraction, discard points outside new hypercube
            % limits
            xlhs = [xlhs(1:nInitialPoints,:); samples];
            gobs = [gobs(1:nInitialPoints,:); gobsSampled];
            fobs = [fobs(1:nInitialPoints,:); fobsSampled];
        end
        
        if norm(xstar(end-1,:) - xstar(end,:))/norm(ub-lb) <= tol2
            
            if largeContractPerformed
                disp('Terminating program, no further improvement possible.');
                disp(xstar(end,:));
                disp(fstar(end));
                
                if ~isempty(fobs(all(gobs <= conTol, 2)))
                    feasibleIndex = find(min(fobs(all(gobs <= conTol, 2))) == fobs);
                    fprintf('\nBest feasible value found: %8.6f\nat point\n',fobs(feasibleIndex))
                    disp(xlhs(feasibleIndex,:));
                    idx = rangesearch(xlhs,xlhs(feasibleIndex,:),0.005*norm(lb-ub));
                    fprintf('%d of points are within 0.5%% euclidian range of this point based on\noriginal domain\n',numel(idx{:})-1);
                    fprintf('Number of funcion evaluations: %d\n',funEvals);
                end
                toc
                break;
            end
            
            if sum(moveTrack(k:(k+1))) >= 2
                % perform large contraction to ensure optimality
                [hub, hlb, exitflagref,ubopt,lbopt] = refineHypercube(xstark,dlb,hlb,dub,hub,0.9999,tolContraction);
                largeContractPerformed = true;
            end
            
        end
        
        xjk = xstark;
        k = k + 1;
        j = 1;
    end
end