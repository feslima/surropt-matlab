function [x, fval, exitflag, output] = minipopt(OBJFUN,CONSFUN,bnds,consbnds,OBJGRAD,CONSJAC,x0)

% sanity check of bounds and constraint bounds
x0 = x0(:)'; % the starting point

[nb_chk, M] = size(bnds);       % M is the number of variables
[nbc_chk, N] = size(consbnds);  % N is the number of constraints

if nb_chk ~= 2 || nbc_chk ~= 2
    error(['You need to specify the variable bounds as a 2-by-M '...
        'matrix where the first line is the lower bound values '...
        'and upper bounds values on the second.']);
end

if M ~= size(x0,2)
    error(['Number of specified bounds of variables is not the '...
        'same as the variables values in the initial estimative.']);
end

% sanity check of the objective and constraint functions handles, also
% their gradients

if ~isa(OBJFUN, 'function_handle') || ~isa(CONSFUN, 'function_handle') || ...
    ~isa(OBJGRAD, 'function_handle') || ~isa(CONSJAC, 'function_handle')
    error (['One or more of the object or constraint (functions and '...
        'gradient/jacobian is not a function handle.']);
end

fval_chk = OBJFUN(x0);

if numel(fval_chk) ~= 1
    error('Objective function must return a scalar value.');
end

cval_chk = CONSFUN(x0);

if numel(cval_chk) ~= N || ~isvector(cval_chk)
    error(['Constraint function evaluation is returning a number of '...
        'different than the number of specified constraint bounds.']);
end

fgrad_chk = OBJGRAD(x0);

if numel(fgrad_chk) ~= M || ~isvector(fgrad_chk)
    error('Objective gradient function must be a vector of M elements');
end

cgrad_chk = CONSJAC(x0);
[ncgrad_chk, mcgrad_chk] = size(cgrad_chk);

if ncgrad_chk ~= N || mcgrad_chk ~= M
    error('The constraint jacobian must be a N-by-M matrix.');
end

% end of sanity check
opts.lb = bnds(1,:);     % Lower bound on the variables
opts.ub = bnds(2,:);     % Upper bound on the variables
opts.cl = consbnds(1,:); % Lower bound on the constraint functions
opts.cu = consbnds(2,:); % Upper bound on the constraint functions

% Initialize the dual point.
% opts.zl     = ones(1,M);
% opts.zu     = ones(1,M);
% opts.lambda = ones(1,N);

% Set the IPOPT options.
% opts.ipopt.mu_strategy = 'adaptive';
opts.ipopt.tol               = 1e-6;
opts.ipopt.constr_viol_tol   = 1e-6;
opts.ipopt.hessian_approximation = 'limited-memory'; % quasi newton hess
% opts.ipopt.print_options_documentation = 'yes';
opts.ipopt.print_level = 0; % turn off report print
% opts.ipopt.option_file_name = 'ipopt.opt';
% opts.ipopt.output_file = 'ipopt.out';

% The callback functions.
funcs.objective         = OBJFUN;
funcs.constraints       = CONSFUN;
funcs.gradient          = OBJGRAD;
funcs.jacobian          = @(x) sparse(CONSJAC(x));
funcs.jacobianstructure = @() sparse(ones(N,M));

% Run IPOPT.
[x, info] = ipopt(x0,funcs,opts);

if nargout >= 2
    fval = OBJFUN(x);
    if nargout >= 3
        exitflag = info.status;
        if nargout >= 4
            output = info;
        end
    end
end

end