function [decMtxA,decVecB] = func_decay(paramUse,nPool,excDt)

T1All = paramUse.T1;
T2All = paramUse.T2;
M0All = paramUse.M0;
mtAll = paramUse.MtRate;

Atemp = zeros(nPool*3,nPool*3);
Btemp = zeros(nPool*3,1);
for ii = 1:nPool
    T2Decay = -1/T2All(ii)-sum(mtAll(ii,:));
    T1Decay = -1/T1All(ii)-sum(mtAll(ii,:));
    Atemp(1+(ii-1)*3,1+(ii-1)*3) = T2Decay;
    Atemp(2+(ii-1)*3,2+(ii-1)*3) = T2Decay;
    Atemp(3+(ii-1)*3,3+(ii-1)*3) = T1Decay;
    Btemp(3+(ii-1)*3) = M0All(ii)/T1All(ii);
    for jj = 1:nPool
        if jj~=ii
            Atemp(1+(jj-1)*3,1+(ii-1)*3) = mtAll(ii,jj);
            Atemp(2+(jj-1)*3,2+(ii-1)*3) = mtAll(ii,jj);
            Atemp(3+(jj-1)*3,3+(ii-1)*3) = mtAll(ii,jj);
        end
    end
end

decMtxA = expm(Atemp*excDt);
decVecB = (decMtxA-eye(nPool*3))*(Atemp\Btemp);


