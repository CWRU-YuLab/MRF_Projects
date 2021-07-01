function [Madc,Mtr] = func_freeprecession(FP,paramUse,Mrf)

nPool = length(Mrf)/3;

%ADC delay
[Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,0,FP.adcDelay);
MpreAdc = Mtx*Mrf+Vec;

%TR
if FP.curDePhase~=0
    paramUseNew = paramUse;
    paramUseNew.OffRes = paramUse.OffRes+FP.curDePhase;
    [Mtx,Vec] = func_rotation_decay(paramUseNew,0,nPool,0,FP.dePhaseDura);
    Mtr = Mtx*MpreAdc+Vec;
else
    [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,0,FP.dePhaseDura);
    Mtr = Mtx*MpreAdc+Vec;
end


%readout
adcDt = FP.adcDura/FP.adcNpt;
if FP.adcNpt>1
    [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,0,adcDt);
    Madc = func_multimult_mex(Mtx,Vec,FP.adcNpt,MpreAdc);
    % Use mex function to calculate following code
    %     for ii = 2:FP.adcNpt
    %         Madc(:,ii) = Mtx*Madc(:,ii-1))+Vec;
    %     end
else
    Madc = MpreAdc;
end

end