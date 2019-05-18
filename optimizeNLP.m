function [x,fval,exitflag] = optimizeNLP(obj_surr, con_surr, x0, lb, ub, solver)

if strcmp(solver, 'ipopt')
    % ipopt as minimizer chosen
    [x,fval,exitflag] = minipopt(@(x) functionObjectivePrediction(x,obj_surr), ...
        @(x) functionConstraintPrediction(x,con_surr),...
        [lb; ub], [-inf(1,numel(con_surr)); zeros(1,numel(con_surr))], ...
        @(x) functionObjectiveGradientPrediction(x,obj_surr),...
        @(x) functionConstraintJacobianPrediction(x,con_surr), x0);
    
    switch exitflag
        case {0, 1}
            exitflag = 1;
        otherwise
            exitflag = -1;; % solver failed
    end
    
elseif strcmp(solver, 'sqp-active-set')
    opts = optimoptions('fmincon','display','none','Algorithm','active-set',...
    'TolCon',conTol,'GradObj','on',...
    'GradConstr','on');
    
    [x,fval,exitflag] = fmincon(@(x) functionObjectivePrediction(x,obj_surr),...
        x0,[],[],[],[],lb,ub,@(x) functionConstraintPrediction(x,con_surr),opts);
else
    error('Invalid NLP solver option.');
end

end

function [f,g] = functionObjectivePrediction(x,krmodelfobj)

[f,g] = predictor(x,krmodelfobj);

end

function g = functionObjectiveGradientPrediction(x,obj_surr)

[~,g] = predictor(x,obj_surr);

end

function [c, ceq, gc, gceq] = functionConstraintPrediction(x,con_surr)

for i = 1:length(con_surr)
    [c(i),gc(:,i)] = predictor(x,con_surr(i));
end

ceq = [];
gceq = [];

end

function jc = functionConstraintJacobianPrediction(x,con_surr)
% returns N-by-M matrix, where N is the number of constraints and M is the
% number of variables
for i = 1:length(con_surr)
    [~,jc(i,:)] = predictor(x,con_surr(i));
end

end