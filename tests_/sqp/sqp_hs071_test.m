function sqp_hs071_test
lb = [1 1 1 1];                            % lower bound
ub = [5 5 5 5];                            % upper bound
x0 = [1 5 5 1];                            % initial estimative
[x, fval, exitflag, nevals, lambda, iter] = sqp(@obj_sqp_test, x0, @nonlcon_sqp_test, lb, ub)
% save('hs071_testvars.mat', 'deb_struct')
% define the constraints function
    function [c,ceq,gc,gceq] = nonlcon_sqp_test(x)
        c = 25 - prod(x);
        ceq = sum(x.^2) - 40;
        gc = - prod(x(:)') ./ x(:)';
        gceq = 2*x(:)';
    end
% define the objective function
    function [f, df] = obj_sqp_test(x)
        f = x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3); % objective function
        df = [ x(1) * x(4) + x(4) * (x(1) + x(2) + x(3)); 
            x(1) * x(4); 
            x(1) * x(4) + 1; 
            x(1) * (x(1) + x(2) + x(3))]; % gradient of objective function
    end
end
