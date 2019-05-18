function octave_sqp_test
x0 = [-1.8 1.7 1.9 -0.8 -0.8];
lb = [];
ub = [];

% [x, fval, exitflag, nevals, lambda, iter, deb_struct] = sqp(@f, x0, @g, lb, ub)

[ci,ce,C,F] = g(x0);
[~,c] = f(x0);
H = eye(5);
[x,fval,exitflag,output,lambda] = quadprog(H,c,C,-ci,F,-ce,[],[])
lambda.eqlin
function [r, req, rgrad, reqgrad] = g(x)
        r = [];
        rgrad = [];
        % equality constraint
        req = [sum(x.^2) - 10;
            x(2) * x(3) - 5*x(4)*x(5);
            x(1)^3 + x(2)^3 + 1];
        reqgrad = [2 .* x;
            0 x(3) x(2) -5*x(5) -5*x(5);
            3*x(1)^2 3*x(2)^2 0 0 0];
    end

    function [obj, objgrad] = f(x)
        obj = exp(prod(x)) - 0.5*(x(1)^3 + x(2)^3 +1)^2;
        objgrad = [prod(x(2:end))*exp(prod(x))-(3*x(1)^2)*(x(1)^3 + x(2)^3 + 1) prod(x(2 ~= 1:5)) * exp(prod(x)) - (3*x(2)^2)*(x(1)^3 + x(2)^3 + 1) prod(x(3 ~= 1:5)) * exp(prod(x)) prod(x(4 ~= 1:5)) * exp(prod(x)) prod(x(5 ~= 1:5)) * exp(prod(x))]';
    end

end