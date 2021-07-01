function [dMtx,dVec] = func_CRLB(paramUse,nPool,excDt)

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

% for the dA/dKfCK
dA_dkf = zeros(nPool*3,nPool*3); 
dA_dkf((1:3),(1:3)) = diag([-1 -1 -1]);
dA_dkf((4:6),(1:3)) = diag([1 1 1]);
dA_dkf((1:3),(4:6)) = diag ([M0All(1)./M0All(2) M0All(1)./M0All(2) M0All(1)./M0All(2)]);
dA_dkf((4:6),(4:6)) = -1.*diag ([M0All(1)./M0All(2) M0All(1)./M0All(2) M0All(1)./M0All(2)]);


AMtx = expm(Atemp*excDt);
dMtx = dA_dkf*expm(Atemp*excDt);
dVec = (AMtx-eye(nPool*3))*(Atemp.^2\Btemp)-(Atemp\Btemp).*AMtx;