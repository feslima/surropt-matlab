% SQP Algorithm - Nocedal/Octave
%
%  Solves NL problem of type:
%                           min f(x)      (objective function)
%                               g(x) <= 0 (NEGATIVE inequality constraints)
%                               h(x) == 0 (equality constraints)
%                           lb <= x <= ub (variables bounds)
%
% INPUTS:
% objfun - function handle that computes the objective and its gradient
% x0 - initial estimative
% con - function handle that computes both equality and inequality
%       constraints. Their jacobian must be included (see example).
% lb - lower bound of independent variables
% ub - upper bound of independent variables
% options - structure that must contain the following fields:
%         - tolopt: optimality tolerance value for the KKT conditions.
%                   Default is 1e-6;
%         - tolstep: minimum step size between two iterations
%                    (step_size = abs(x_new - x_old). Default is 1e-6;
%         - tolcon: maximum constraint violation. Default is 1e-6;
%         - maxiter: maximum number of iterations. Default is 400;
%         - maxfunevals: maximum number of function evaluations.
%                        Default is 100 * number of elements in the initial
%                        estimative;
%
% OUTPUTS:
% x - optimum value found if exitflag > 0 (i.e. problem has converged)
% obj - function objective value found at x
% exiflag - 1 (optimum was found);
%           2 (step size is too small and constraint violation is less than
%              min constraint violation.);
%           0 (too many iterations - did not converge);
%         - < 0 (negative - infeasible. NEEDS DOCUMENTATION FINISH);
%
% EXAMPLE:
% % define the constraints function
% function [c,ceq,gc,gceq] = nonlcon_sqp_test(x)
% c = 25 - prod(x);
% ceq = sum(x.^2) - 40;
%
% gc = - prod(x(:)') ./ x(:)';
% gceq = 2*x(:)';
% end
%
% % define the objective function
% function [f, df] = obj_sqp_test(x)
%     f = x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3); % objective function
%     df = [ x(1) * x(4) + x(4) * (x(1) + x(2) + x(3));
%         x(1) * x(4);
%         x(1) * x(4) + 1;
%         x(1) * (x(1) + x(2) + x(3))]; % gradient of objective function
% end
%
% lb = [1 1 1 1]                            % lower bound
% ub = [5 5 5 5]                            % upper bound
% x0 = [1 5 5 1]                            % initial estimative
%
% [x, fval, exitflag, nevals, lambda, iter] = sqp(@obj_sqp_test, x0, @nonlcon_sqp_test, lb, ub)

function [x, obj, exitflag, nevals, lambda, iter, deb_struct] = sqp(...
    objfun, x0, confun, lb, ub, options)
% TODO: Error treatment.; (11/02/2019) DONE: (14/02/2019)
% NEEDS: - function handle check (see if they are proper handles that
%          enable calculations);
%        - verify if x0, lb and ub are vector of same size;
%        - verify if x0 is between lb and ub;
%        - verify if there are inconsistencies with lb and ub
%          (e.g. lb > ub); DONE: (13/02/2019)

% TODO: error/routine check when one of either equality or inequality comes
% as an empty handle. (Proposal: make it a empty function when this
% happens); (12/02/2019) DONE (12/02/2019);

% TODO: transform bounds into inequality constraints.; (12/02/2019)
% DONE (12/02/2019)

% TODO: implement maxfunevals check and exitflag.; (12/02/2019)
% DONE: (13/02/2019)

% Check input vectors
if isvector(x0) && isvector(lb) && isvector(ub)
    x0 = x0(:); % force to be column
    lb = lb(:);
    ub = ub(:);
    
    if numel(x0) ~= numel(lb) || numel(x0) ~= numel(ub)
        error(['Initial estimative (x0) and bounds (lb/ub) '...
            'do not have the same size.']);
    end
    
    if any(lb >= ub) % check if lb >= ub
        error('Lower bound must be lower than upper bound (lb < ub).');
    end
else
    error('x0, lb and ub must be vectors');
end

if nargin > 9 % options structure is set
    if ~isstruct(options) % if options is not a struct, throw exception
        error('options must be a struct variable.');
    else
        % check if the options structure has the required fields
        
        % maximum number of iterations
        if isfield(options,'maxiter')
            if isscalar(options.maxiter) && options.maxiter > 0 ...
                    && fix(options.maxiter) == options.maxiter
                % if is a positive scalar integer, set the value
                maxiter = options.maxiter;
            else
                error(...
                    'Invalid value set for maximum number of iterations.');
            end
        else
            maxiter = 400; % set to default
        end
        
        % maximum number of function evaluations
        if isfield(options,'maxfunevals')
            if isscalar(options.maxfunevals) && options.maxfunevals > 0 ...
                    && fix(options.maxfunevals) == options.maxfunevals
                % positive scalar integer
                maxfunevals = options.maxfunevals;
            else
                error(...
                    ['Invalid value set for maximum number ' ...
                    'of function evaluations.']);
            end
        else
            maxfunevals = 100 * numel(x0);
        end
        
        % minimum step size
        if isfield(options,'tolstep')
            if isscalar(options.tolstep) && options.tolstep > 0
                % if is a positive scalar, set it
                tolstep = options.tolstep;
            else
                error('Invalid value set for step size tolerance.');
            end
        else
            tolstep = 1e-6;
        end
        
        % optimality tolerance
        if isfield(options,'tolopt')
            if isscalar(options.tolopt) && options.tolopt > 0
                % if is a positive scalar, set it
                tolopt = options.tolopt;
            else
                error('Invalid value set for tolopt.');
            end
        else
            tolopt = 1e-6;
        end
        
        % constraint tolerance
        if isfield(options,'tolcon')
            if isscalar(options.tolcon) && options.tolcon > 0
                % if it is a positive scalar, set it
                tolcon = options.tolcon;
            else
                error('Invalid value set for constraint tolerance.');
            end
        else
            tolcon = 1e-6;
        end
    end
    
else % create the default parameters
    maxiter = 400;
    maxfunevals = 100 * numel(x0);
    tolstep = 1e-6;
    tolopt = 1e-6;
    tolcon = 1e-6;
end

% check if function handlea are working properly:
% - objective function and gradient must return scalar and vector
if isa(objfun,'function_handle')
    
    % if is a function handle
    [ofun_eval, ogfun_eval] = objfun(x0);
    
    if ~isscalar(ofun_eval) || ~isnumeric(ofun_eval)
        % if is not scalar, empty or not numeric throw error
        error(['Objective function must return a valid '...
            'numeric scalar value.']);
    end
    
    if ~iscolumn(ogfun_eval) || ~isnumeric(ogfun_eval)
        % if it is not a vector, empty or not numeric throw error
        error(['Objective function gradient must return '...
            'a valid numeric column vector.']);
    end
    
else
    error('Objective function and gradient must be a function handle.');
end

% check if any of the constraints are coming as an empty handle and make
% them an empty function
if isempty(confun)
    confun = @(x) empty_nonlcon(x); % make it an empty handle
    
elseif isa(confun,'function_handle')
    % check constraint and their jacobian
    [cif_eval, cef_eval, cigf_eval, cegf_eval] = confun(x0);
    
    % check if constraint functions are returning column vectors
    if ~iscolumn(cif_eval) && ~isempty(cif_eval)
        error(['Inequality constraint function must return a column '...
            'vector.']);
    end
    
    if ~iscolumn(cef_eval) && ~isempty(cef_eval)
        error(['Equality constraint function must return a column '...
            'vector.']);
    end
    
    % check if the jacobian is specified without its corresponding handle
    if (isempty(cif_eval) && ~isempty(cigf_eval)) || ...
            (isempty(cef_eval) && ~isempty(cegf_eval))
        error(['You can''t specify jacobian without their '...
            'constraint functions.']);
    end
    
    % check if the jacobians have the proper dimensions
    [n_ce, n_var] = size(cegf_eval);
    if n_ce ~= 0 && (n_ce ~= numel(cef_eval) || n_var ~= numel(x0))
        error(['The jacobian of equality constraints must have '...
            'its number of rows equal to the number of equality '...
            'constraints functions, and its number of columns equal '...
            'to the number of elements of the initial estimate.']);
    end
    
    [n_ci, n_var] = size(cigf_eval);
    if n_ci ~= 0 && (n_ci ~= numel(cif_eval) || n_var ~= numel(x0))
        error(['The jacobian of inequality constraints must have '...
            'its number of rows equal to the number of inequality '...
            'constraints functions, and its number of columns equal '...
            'to the number of elements of the initial estimate.']);
    end
else
    error(['Constraints function and their jacobian must '...
        'be a function handle.']);
end

% Make bound into constraints
% lbidx = true(size(x0)); ubidx = lbidx;  % instantiation (necessary??)
lbgrad = - eye(numel(x0)); ubgrad = -lbgrad;

lbidx = lb ~= -Inf; % get index of bounds different than +/- infinity
ubidx = ub ~=  Inf;
lbgrad = lbgrad(lbidx, :);
ubgrad = ubgrad(ubidx, :);

confun = @(x) bnd2cf(x, lbidx, ubidx, lb, ub, [lbgrad; ubgrad], confun);
% cigf = @(x) bnd2cgf(x, [lbgrad; ubgrad], cigf);

% global structure for parametrization
globals = struct;

% Initialize variables: objective and constraints (gradients and jacobians)
x = x0;

[obj, c] = objfun(x0);
% obj = ofun(x0);
% c = ogfun(x0);
globals.nevals = 1;

% Initialize the positive definite Hessian approximation
n = numel(x0);
B = eye(n);

% Evaluate the constraint info
[ci, ce, C, F] = confun(x0);
% ce = cef(x0);
% F = cegf(x0);
%
% ci = cif(x0);
% C = cigf(x0);

A = [F; C]; % Jacobian matrix of both equalities and inequalities

% Choose an initial lambda

lambda = 100 * ones(numel(ce) + numel(ci), 1);

iter = 0;

% turn off warning for not symmetric hessian in the QP solver
warning('off','optim:quadprog:HessianNotSym');

p = [];
while iter < maxiter
    % Convergence check
    nr_f = numel(ce); % number of equality constraints
    
    % split the lagrange multipliers
    %     lambda_e = lambda(1:nr_f);
    lambda_i = lambda((nr_f+1):end);
    
    con = [ce; ci];
    
    % Karush-Kuhn-Tucker conditions
    % TODO: Implement MATLAB's infinity norm check to improve check
    t0 = norm(c - A' * lambda);
    t1 = norm(ce);
    t2 = all(ci <= 0);
    t3 = all(lambda_i >= 0);
    t4 = norm(lambda .* con);
    
    if (t2 && t3 && max([t0 t1 t4]) < tolopt)
        % Problem has converged with all constraints being satisfied.
        exitflag = 1;
        break
    end
    
    % Solve the QP subproblem to compute the search direction p
    lambda_old = lambda; % store old multipliers
    pold = p;
    [p, ~, qpexflag, ~, lmbda_struct] = quadprog(B, c, C, -ci, F, -ce, ...
        [],[],x, optimoptions('quadprog','Display','off'));
    
%     deb_struct{iter+1}.x = x;
%     deb_struct{iter+1}.B = B;
%     deb_struct{iter+1}.p = p;
%     deb_struct{iter+1}.exflag = qpexflag;
%     deb_struct{iter+1}.lmbda = [lmbda_struct.eqlin; lmbda_struct.ineqlin];
    % for constraints of type g(x) <= 0 use the syntax below
    %     [p, ~, qpexflag, ~, lmbda_struct] = quadprog(B, c, C, -ci, F, -ce, ...
    %         [],[],[], optimoptions('quadprog','Display','off'));
    
    % QP solution check
    if qpexflag == 0
        warning(['QP failed to find a solution'...
            '(Maximum number of iterations reached).']);
    elseif qpexflag == -3
        warning('QP is unbounded.');
    elseif qpexflag < 0
        warning('QP is infeasible')
        lambda = lambda_old;
        p = x;
        
    elseif qpexflag > 0 % no problem
        % update multipliers
        lambda(1:nr_f) = -lmbda_struct.eqlin;
        lambda((nr_f+1):end) = -lmbda_struct.ineqlin;
    end
    
    % perform linesearch
    [x_new, ~, obj_new, c_new, globals] = ...
        linesearch(x, p, objfun, confun, lambda, obj, c, globals);
    
    % Re-evaluate objective, constraints and gradients at new x value.
    [ci_new, ce_new, C_new, F_new] = confun(x_new);
    
    %     ce_new = cef(x_new);
    %     F_new = cegf(x_new);
    %
    %     ci_new = cif(x_new);
    %     C_new = cigf(x_new);
    
    A_new = [F_new; C_new];
    
    % set s and y for BFGS update
    y = c_new - c;
    
    if ~isempty(A)
        t = (A_new - A)' * lambda;
        y = y - t;
    end
    
    delx = x_new - x; % step size
    
    % Check if step size is too small
    if norm(delx) < tolstep * norm(x)
        % check for minimum constraint violation
        if max([ce_new; ci_new;]) < tolcon
            exitflag = 2; % Step size is too small and constraint violation
            % is less than min constraint violation.
            
        else
            exitflag = -2; % Step size is too small but there is some
            % significant constraint violation.
            
        end
        break;
    end
    
    % check for number of function evaluations
    if globals.nevals > maxfunevals
        exitflag = 0;
        break;
    end
    
    % perform the damped BFGS update
    delxt = delx';
    
    d1 = delxt * B * delx;
    
    t1 = 0.2 * d1;
    t2 = delxt * y;
    
    if t2 < t1
        theta = 0.8 * d1 / (d1 - t2);
    else
        theta = 1;
    end
    
    r = theta * y + (1 - theta) * B * delx;
    
    d2 = delxt * r;
    
    % if d1 or d2 -> 0, the BFGS update will fail.
    if (d1 == 0 || d2 == 0)
        exitflag = -1; % BFGS update failed
        break
    end
    
    B = B - B*delx*delxt*B / d1 + r*r' / d2;
    
    % update values
    x = x_new;
    obj = obj_new;
    c = c_new;
    
    ce = ce_new;
    F = F_new;
    
    ci = ci_new;
    C = C_new;
    
    A = A_new;
    
    iter = iter + 1;
end % while iter < maxiter

if iter >= maxiter
    exitflag = 0; % too many iterations
end

nevals = globals.nevals;
x = x'; % force to be row vector

end % sqp

% Linesearch auxiliary routine
function [x_new, alpha, obj, objgrad, globals] = ...
    linesearch(x, p, objfun, nonlcon, lambda, obj, objgrad, globals)
% choose eta and tau ([0, 0.5] and [0, 1] are their respective ranges).
eta = 0.25;
tau = 0.5;

% Choose mu to satisfy 18.36 with sigma = 1
delta_bar = sqrt(eps);

if isempty(lambda) % no constraint multipliers
    mu = 1 / delta_bar;
else
    mu = 1 / (norm(lambda, Inf) + delta_bar);
end

alpha = 1;

c = objgrad;
% [~, c] = objfun(x);

[d, ce] = nonlcon(x);
% d = cif(x);
% ce = cef(x);

% call the merit function
[phi_x_mu, obj, objgrad, globals] = merit_fun(obj, objgrad, objfun, ...
    nonlcon, x, mu, globals);

D_phi_x_mu = c' * p; % directional derivative

% only those elements of d corresponding to violated constraints should be
% included.
idx = d > 0;
t = - norm([ce; d(idx)],1) / mu;
if ~isempty(t)
    D_phi_x_mu = D_phi_x_mu + t;
end

while true % while loop to find alpha that causes decrease in merit
    [p1, obj, objgrad, globals] = merit_fun([], objgrad, objfun, ...
        nonlcon, x + alpha*p, mu, globals);
    p2 = phi_x_mu + eta * alpha * D_phi_x_mu;
    
    if p1 > p2
        % Reset alpha = tau_alpha * alpha for some tau_alpha in the range
        % [0, 1].
        tau_alpha = 0.9 * tau; % why 0.9? is it a general value?
        alpha = tau_alpha * alpha;
    else
        break % the alpha is enough to cause decrease
    end
end % while true

x_new = x + alpha * p;

end % linesearch

% Merit function evaluation auxiliary routine
function [merit, obj, objgrad, globals] = merit_fun(obj, objgrad, ...
    objfun, nonlcon, x, mu, globals)

[ci, ce] = nonlcon(x);
% ce = cef(x);
% ci = cif(x);

idx = ci > 0; % indexes of negative (violated) constraints

con = [ce; ci(idx)];

if isempty(obj)
    [obj, objgrad] = objfun(x); % evaluate function if obj value is empty
    globals.nevals = globals.nevals + 1;
end

merit = obj;
t = norm(con,1) / mu;

if ~isempty(t) % if t is not empty
    merit = merit + t;
end

end % merit_fun

% Transform bounds into constraints (gradients included)
function [cif, cef, cigf, cegf] = bnd2cf(x, lbidx, ubidx, lb, ub, ...
    bnds_grad, nonlcon)

[cif_eval, cef, cigf_eval, cegf] = nonlcon(x);
if ~isempty(cif_eval)
    cif = [cif_eval; lb(lbidx) - x(lbidx); x(lbidx) - ub(ubidx)];
else
    cif = [lb(lbidx) - x(lbidx); x(lbidx) - ub(ubidx)];
end

if ~isempty(cigf_eval)
    cigf = [cigf_eval; bnds_grad];
else
    cigf = bnds_grad;
end

end % bnd2cf

% Empty non-linear constraint function (if nonlcon is empty)
function [cif, ceq, gcif, gceq] = empty_nonlcon(~)
cif = []; ceq = cif; gcif = cif; gceq = cif;
end % empty_nonlcon