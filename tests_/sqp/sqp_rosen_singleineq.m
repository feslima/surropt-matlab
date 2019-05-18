function sqp_rosen_singleineq
x0 = [-1 2];

[ci,ce,C,F] = rosen_1con(x0);
[~,c] = rosenbrockwithgrad(x0);
H = eye(numel(x0));

[x,fval,exitflag,output,lambda] = quadprog(H,c,C,-ci,F,-ce,[],[])

function [f,g] = rosenbrockwithgrad(x)
        % Calculate objective f
        f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
        
        if nargout > 1 % gradient required
            g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
                200*(x(2)-x(1)^2)];
        end
    end

    function [c, ceq, cgrad, ceqgrad] = rosen_1con(x)
        c = x(1) + 2*x(2) - 1;
        cgrad = [1 2];
        ceq = [];
        ceqgrad = [];
    end
end