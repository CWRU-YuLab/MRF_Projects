function [Mtx,Vec] = func_rotation_decay(paramUse,b1,nPool,gamma,dt)

T1All = paramUse.T1;
T2All = paramUse.T2;
M0All = paramUse.M0;
offResAll = -paramUse.OffRes;
mtAll = paramUse.MtRate;

vecTemp = zeros(nPool*3,1);
vecTemp(3:3:end) = M0All./T1All;

mtxTemp = zeros(nPool*3,nPool*3);
T1Decay = reshape(-1./T1All,[],1)-sum(mtAll,2);
T2Decay = reshape(-1./T2All,[],1)-sum(mtAll,2);
for rid = 1:nPool
    for cid = 1:nPool
        if(rid~=cid)
            mtxTemp((1:3)+(rid-1)*3,(1:3)+(cid-1)*3) = diag(repmat(mtAll(cid,rid),[1 3]));
        else
            mtxTemp((1:3)+(rid-1)*3,(1:3)+(rid-1)*3) = diag([T2Decay(rid) T2Decay(rid) T1Decay(rid)]);
            mtxTemp(1+(rid-1)*3,2+(rid-1)*3) = -2*pi*offResAll(rid);
            mtxTemp(2+(rid-1)*3,1+(rid-1)*3) = -mtxTemp(1+(rid-1)*3,2+(rid-1)*3);
            mtxTemp(1+(rid-1)*3,3+(rid-1)*3) = -2*pi*gamma*imag(b1);
            mtxTemp(3+(rid-1)*3,1+(rid-1)*3) = -mtxTemp(1+(rid-1)*3,3+(rid-1)*3);
            mtxTemp(2+(rid-1)*3,3+(rid-1)*3) = 2*pi*gamma*real(b1);
            mtxTemp(3+(rid-1)*3,2+(rid-1)*3) = -mtxTemp(2+(rid-1)*3,3+(rid-1)*3);
        end
    end
end

Mtx = expm(mtxTemp*dt);
Vec = (Mtx-eye(nPool*3))*(mtxTemp\vecTemp);
