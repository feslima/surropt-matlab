% hlb - hypercube lower bound
% hub - hypercube upper bound
% dlb - domain lower bound
% dub - domain upper bound

% output:
% exitflag = 1 (contract inside the domain - figure 4.a)
% exitflag = 2 (move inside the domain - figure 4.b)
% exitflag = 3 (contract at the domain limit - figure 4.c)
% exitflag = 4 (move at the domain limit - figure 4.d)

function [newUB, newLB, exitflag,UBopt,LBopt] = refineHypercube(xstark,dlb,hlb,dub,hub,contractFactor,tolContraction)

% Is (xstar)k inside a limit of sampling hypercube
range = hub - hlb; % hypercube range

% Check if xstark is at domain bound
if isInsideHyperCube(xstark,dlb,dub) == 1 % no
    
    if isInsideHyperCube(xstark,hlb,hub) == 1 && norm(hlb-hub)/norm(dlb-dub) >= tolContraction %inside hypercube
        % center and contract
        newUB = xstark + (1 - contractFactor).*range./2;
        newLB = xstark - (1 - contractFactor).*range./2; 
        exitflag = 1;
    elseif isInsideHyperCube(xstark,hlb,hub) ~= 1 || norm(hlb-hub)/norm(dlb-dub) < tolContraction %at limit
        % center and move
        newUB = xstark + range./2;
        newLB = xstark - range./2;
        exitflag = 2;
    end
    
    LBopt = newLB;
    UBopt = newUB;
    if any(bsxfun(@lt,newLB,dlb))        
        LBopt(bsxfun(@lt,newLB,dlb)) = dlb(bsxfun(@lt,newLB,dlb));
    end
    
    if any(bsxfun(@gt,newUB,dub))        
        UBopt(bsxfun(@gt,newUB,dub)) = dub(bsxfun(@gt,newUB,dub));
    end
    
elseif isInsideHyperCube(xstark,dlb,dub) ~= 1 % no (at domain limit)
    
    % Check if xstark is at hypercube bound
    if isInsideHyperCube(xstark,hlb,hub) == 1 && norm(hlb-hub)/norm(dlb-dub) >= tolContraction %inside hypercube
        % center and contract
        newUB = xstark + (1 - contractFactor).*range./2;
        newLB = xstark - (1 - contractFactor).*range./2;
        exitflag = 3;
    elseif isInsideHyperCube(xstark,dlb,dub) ~= 1 || norm(hlb-hub)/norm(dlb-dub) < tolContraction %at hypercube limit
        % center and move
        newUB = xstark + range./2;
        newLB = xstark - range./2;
        exitflag = 4;
    end
    
    LBopt = newLB;
    UBopt = newUB;
    if any(bsxfun(@lt,newLB,dlb))        
        LBopt(bsxfun(@lt,newLB,dlb)) = dlb(bsxfun(@lt,newLB,dlb));
    end
    
    if any(bsxfun(@gt,newUB,dub))        
        UBopt(bsxfun(@gt,newUB,dub)) = dub(bsxfun(@gt,newUB,dub));
    end
end