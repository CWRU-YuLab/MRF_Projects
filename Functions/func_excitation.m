function [Mout Mout2] = func_excitation(RF,paramUse,b1Scale,gamma,M)

%prepare RF shape stored in excB1, the timing is dependent on excDt and
%RF.Npt
if strcmp(RF.Name,'GaussTrunc20')
    %excWin = reshape(gausswin(RF.Npt,1.8),[],1); %a bit time consuming
    if(mod(RF.Npt,2)==0)
        RF.Npt = RF.Npt+1;
    end
    alpha = 1.8;
    n = -(RF.Npt-1)/2:(RF.Npt-1)/2;
    stdev = (RF.Npt-1)/2/alpha;
    excWin = exp(-0.5*(n/stdev).^2);
    excDt = RF.Dura/RF.Npt;
    excB1max = RF.FA*b1Scale/180*pi/(sum(excWin*excDt)*gamma*2*pi);
    excB1 = excB1max*excWin*exp(1i*RF.Ph/180*pi); %with RF phase applied
elseif strcmp(RF.Name,'bp')
    excWin = ones(RF.Npt,1);
    excDt = RF.Dura/RF.Npt;
    excB1max = RF.FA*b1Scale/180*pi/(sum(excWin*excDt)*gamma*2*pi);
    excB1 = excB1max*excWin*exp(1i*RF.Ph/180*pi); %with RF phase applied
elseif strcmp(RF.Name,'sinc7h')
    RF_TBW = 14;
    excDt = RF.Dura/RF.Npt;
    tt = -RF.Dura/2:excDt:(RF_Dur/2-excDt);
    excWin = (1+cos(2*pi.*tt/RF_Dur))/2.*sinc(RF_TBW.*tt/RF_Dur);
    excB1max = RF.FA*b1Scale/180*pi/(sum(excWin*excDt)*gamma*2*pi);
    excB1 = excB1max*excWin*exp(1i*RF.Ph/180*pi);
elseif strcmp(RF.Name,'sech')
    load('InvSech.mat');
    fmag = fmag./max(fmag);
    excWin = fmag(1:8:end).*exp(1i*fph(1:8:end)/180*pi);
    excDt = RF.Dura/length(excWin);
    RF.Npt = length(excWin);
    excB1max = RF.FA*b1Scale/180*pi/(gamma*2*pi*RF.Dura)/(3.162277e-02);
    excB1 = excB1max*excWin*exp(1i*RF.Ph/180*pi);
end

nPool = length(M)/3;

if strcmp(RF.Name,'bp')
    paramUseNew = paramUse;
    paramUseNew.OffRes = paramUse.OffRes-RF.OffRes;
    [Mtx,Vec] = func_rotation_decay(paramUseNew,excB1(1),nPool,gamma,RF.Dura);
    Mout = Mtx*M+Vec;
else
    [decMtxA,decMtxB] = func_decay(paramUse,nPool,excDt); %relaxation and MT effect
    [dMtx,dVec] = func_CRLB(paramUse,nPool,excDt); %relaxation and MT effect
    rX = 2*pi*gamma*real(reshape(excB1,1,[]))*excDt;
    rY = 2*pi*gamma*imag(reshape(excB1,1,[]))*excDt;
    rZ = reshape(2*pi*(RF.OffRes-paramUse.OffRes)*excDt,[],1);
    MoutTemp = func_excitation_mex(rX,rY,rZ,decMtxA,decMtxB,M);
    MoutTemp2 = func_excitation_mex(rX,rY,rZ,dMtx,dVec,M); % @@ added 
    Mout = MoutTemp(:,end);
    Mout2 = MoutTemp2(:,end);
    % The mex function did the following calculation
    %    for ii = 1:RF.Npt
    %         rotMtx = func_rotation(paramUse.OffRes,excB1(ii),nPool,gamma,excDt);
    %         Mout = decMtxA*(rotMtx*Mout)+decMtxB;      
    %    end    
end

%correct the phase relative to reference rotating frame
if(RF.OffRes~=0)
    rotMtx = zeros(nPool*3,nPool*3);
    rotAng = 2*pi*RF.OffRes*RF.Dura;
    for ii = 1:nPool
        rotMtx((1:2)+(ii-1)*3,(1:2)+(ii-1)*3) = [cos(rotAng) sin(rotAng); -sin(rotAng) cos(rotAng)];
        rotMtx(3+(ii-1)*3,3+(ii-1)*3) = 1;
    end
    Mout = rotMtx*Mout;
    %Mout2 = rotMtx*Mout2;
end
end