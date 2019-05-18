function [krmodelfobj,krmodelcon,krmodel,perf] = buildSurrogate(designData,surrObsData,...
    regModel,corrModel,krmodelfObjPrev,krmodelConPrev)

nVDes = size(designData,2);
one = ones(1,nVDes);
theta0 = one;
lbtheta = 1e-12*one;
ubtheta = 100*one;

nSur = size(surrObsData,2);

if nargin <= 4 % if krmodel NOT specified, optimize hyperparameters
    for i = 1:nSur
        [krmodel(i),perf(i)] = dacefit(designData,surrObsData(:,i),regModel,corrModel,...
            theta0,lbtheta,ubtheta);
    end
else % krmodel given, do NOT optimize hyperparameters
    krmodelPrev(1) = krmodelfObjPrev;
        for i = 2:(length(krmodelConPrev)+1)
            krmodelPrev(i) = krmodelConPrev(i-1);
        end
    for i = 1:nSur
        [krmodel(i),perf(i)] = dacefit(designData,surrObsData(:,i),regModel,corrModel,...
            krmodelPrev(i).theta);
    end
end


krmodelfobj = krmodel(1);
krmodelcon = krmodel(2:end);