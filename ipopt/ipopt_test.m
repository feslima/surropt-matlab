% minipopt test function (h071 formulation)
function ipopt_test

x0 = [ 1 5 5 1];
bnds = [1 1 1 1; 5 5 5 5];
cbnds = [25 40; inf 40];

objfun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
consfun = @(x) [ prod(x); sum(x.^2) ];
consjac = @(x) sparse([ prod(x)./x; 2*x ]);

[x, fval, exitflag, output] = minipopt(objfun,consfun,bnds,cbnds,@objgrad,consjac,x0)

% hs051

% ----------------------------------------------------------------------
    function g = objgrad(x)
        g = [ x(1)*x(4) + x(4)*sum(x(1:3))
            x(1)*x(4)
            x(1)*x(4) + 1
            x(1)*sum(x(1:3)) ];
    end
end