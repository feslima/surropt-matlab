% evaporator model that acceps 4 or 7 values as input vector. If the vector
% has size 7, the first 3 elements are the disturbances, the last 4 are the
% MV's. If the vector has 4 elements, all of the are just MV's.

function doeBuild = evaporatorDOE(x)

if size(x,2) == 7
    casos = x;      % X1 T1 T200 F1 F3 P100 F200
    no = x(:,1:3);
elseif size(x,2) == 4
    casos = x;      % F1 F3 P100 F200
    no = [5 40 25]; % X1   T1 T200
else
    error('x must have 4 or 7 columns.');
end

options = optimoptions('fsolve','Display','none');
for i = 1:size(casos,1)
    if size(x,2) > 4
        X1_LHS = casos(i,1); T1_LHS = casos(i,2); T200_LHS = casos(i,3);
        F1_LHS = casos(i,4); F3_LHS   = casos(i,5);
        P100_LHS = casos(i,6); F200_LHS = casos(i,7);
        no = casos(i,1:3);
    else
        X1_LHS = no(1); T1_LHS = no(2); T200_LHS = no(3);
        F1_LHS = casos(i,1); F3_LHS   = casos(i,2);
        P100_LHS = casos(i,3); F200_LHS = casos(i,4);
    end
    [out,fval,exitflag,output] = fsolve(@NonOptimizationModel,1000*ones(12,1),options);
    endIndex = i;
    %     check_matrix(endIndex + 1,1) = exitflag; check_matrix(endIndex + 1,2)= sum(fval.^2);
    doeBuild(endIndex,1:3) = no;
    if size(x,2) > 4
        doeBuild(endIndex,1:7) = casos(i,:);
    else
        doeBuild(endIndex,4:7) = casos(i,:);
    end
    doeBuild(endIndex,8:19) = out(:)';
    doeBuild(endIndex,20) = 600*out(8) + 0.6*F200_LHS + 1.009*(out(1) + F3_LHS) + 0.2*F1_LHS - 4800*out(1);
end

% save([doepath matfilepath],'doeBuild')
    function feq = NonOptimizationModel(x)
        F1 = F1_LHS; F2 = x(1); F3 = F3_LHS; F4 = x(2); F5 = x(3);
        X1 = X1_LHS; X2 = x(4); T1 = T1_LHS; T2 = x(5); T3 = x(6);
        %         L2 = x(7);
        P2 = x(7); F100 = x(8); T100 = x(9); P100 = P100_LHS;
        Q100 = x(10); F200 = F200_LHS; T200 = T200_LHS; T201 = x(11); Q200 = x(12);
        
        feq(1) = (F1 - F4 - F2)/20;
        feq(2) = (F1 * X1 - F2 * X2)/20;    % 39
        feq(3) = (F4 - F5)/ 4;              % 40
        
        feq(4) = (0.5616 * P2 + 0.3126 * X2 + 48.43 - T2);    % 41
        feq(5) = (0.507 * P2 + 55 - T3);                      % 42
        feq(6) = ((Q100 - 0.07 * F1 * (T2 - T1))/38.5 - F4);  % 43
        feq(7) = (0.1538 * P100 + 90 - T100);                 % 44
        feq(8) = (0.16 * (F1 + F3) * (T100 - T2) - Q100);     % 45
        feq(9) = (Q100/ 36.6 - F100);                         % 46
        quoc = (T3 - T200) /( 0.14 * F200 + 6.84 );
        feq(10) = (0.9576 * F200 * quoc - Q200);              % 47
        feq(11) = (T200 + 13.68 * quoc - T201);               % 48
        feq(12) = (Q200/38.5 - F5);                           % 49
        
    end

end