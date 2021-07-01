function Mprep = func_preparation(PreparationInfo,paramUse,b1Scale,gamma,M)

nPrepmodule = length(PreparationInfo);
nPool = length(M)/3;
Mprep = M;
for ii = 1:nPrepmodule
    if strcmp(PreparationInfo{ii}.Type,'InversionSech')
        RF.Name = 'sech';
        RF.FA = PreparationInfo{ii}.InvFA;
        RF.Ph = PreparationInfo{ii}.InvPh;
        RF.OffRes = PreparationInfo{ii}.InvOffRes;
        RF.Dura = PreparationInfo{ii}.InvDura;
        RF.Npt = 256; %Not used here
        Mprep = func_excitation(RF,paramUse,b1Scale,gamma,Mprep);
        [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,gamma,PreparationInfo{ii}.InvDelay);
        Mprep = Mtx*Mprep+Vec;
    elseif strcmp(PreparationInfo{ii}.Type,'Inversion')
        RF.Name = PreparationInfo{ii}.InvName;
        RF.FA = PreparationInfo{ii}.InvFA;
        RF.Ph = PreparationInfo{ii}.InvPh;
        RF.OffRes = PreparationInfo{ii}.InvOffRes;
        RF.Dura = PreparationInfo{ii}.InvDura;
        RF.Npt = PreparationInfo{ii}.InvNpt; %Not used here
        Mprep = func_excitation(RF,paramUse,b1Scale,gamma,Mprep);
        [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,gamma,PreparationInfo{ii}.InvDelay);
        Mprep = Mtx*Mprep+Vec;
    elseif strcmp(PreparationInfo{ii}.Type,'Bp-Saturation')
        RF.Name = 'bp';
        RF.FA = PreparationInfo{ii}.SatB1*2*pi*gamma*PreparationInfo{ii}.SatDura/pi*180;
        RF.Ph = PreparationInfo{ii}.SatPh;
        RF.OffRes = PreparationInfo{ii}.SatOffRes;
        RF.Dura = PreparationInfo{ii}.SatDura;
        RF.Npt = 51;
        Mprep = func_excitation(RF,paramUse,b1Scale,gamma,Mprep);
    elseif strcmp(PreparationInfo{ii}.Type,'Crusher')
        [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,gamma,PreparationInfo{ii}.CrushDura);
        Mprep = Mtx*Mprep+Vec;
        Mprep(1:3:end) = 0;
        Mprep(2:3:end) = 0;
    elseif strcmp(PreparationInfo{ii}.Type,'Delay')
        [Mtx,Vec] = func_rotation_decay(paramUse,0,nPool,gamma,PreparationInfo{ii}.DelayDura);
        Mprep = Mtx*Mprep+Vec;
    elseif strcmp(PreparationInfo{ii}.Type,'T2-MLEV4')
        RF90.Name = 'bp';
        RF90.FA = 90;
        RF90.Ph = PreparationInfo{ii}.ExcPh;
        RF90.OffRes = PreparationInfo{ii}.ExcOffRes;
        RF90.Dura = PreparationInfo{ii}.ExcDura;
        RF90.Npt = 6;
        RF240 = RF90;
        RF240.FA = 240;
        RF240.Ph = RF90.Ph+90;
        RF240.Dura = RF90.Dura*240/90;
        RF240.Npt = 16;
        RF90n = RF90;
        RF90n.Ph = RF90.Ph+180;
        RF240n = RF240;
        RF240n.Ph = RF240.Ph+180;
        TM8 = PreparationInfo{ii}.MixTime/8;
        TM4 = PreparationInfo{ii}.MixTime/4;
        [MtxTM8,VecTM8] = func_rotation_decay(paramUse,0,nPool,gamma,TM8);
        [MtxTM4,VecTM4] = func_rotation_decay(paramUse,0,nPool,gamma,TM4);
        Mprep = func_excitation(RF90,paramUse,b1Scale,gamma,Mprep); %90x
        Mprep = MtxTM8*Mprep+VecTM8;
        for jj = 1:2
            Mprep = func_excitation(RF90,paramUse,b1Scale,gamma,Mprep); %90x
            Mprep = func_excitation(RF240,paramUse,b1Scale,gamma,Mprep); %240y
            Mprep = func_excitation(RF90,paramUse,b1Scale,gamma,Mprep); %90x
            Mprep = MtxTM4*Mprep+VecTM4;
        end
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep);
        Mprep = func_excitation(RF240n,paramUse,b1Scale,gamma,Mprep);
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep);
        Mprep = MtxTM4*Mprep+VecTM4;
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep);
        Mprep = func_excitation(RF240n,paramUse,b1Scale,gamma,Mprep);
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep);
        Mprep = MtxTM8*Mprep+VecTM8;
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep); %90-x
    elseif strcmp(PreparationInfo{ii}.Type,'T2-Adiabatic')
        RF90.Name = 'bp';
        RF90.FA = 90;
        RF90.Ph = PreparationInfo{ii}.ExcPh;
        RF90.OffRes = PreparationInfo{ii}.ExcOffRes;
        RF90.Dura = PreparationInfo{ii}.ExcDura;
        RF90.Npt = 6;
        RF90n = RF90;
        RF90n.Ph = RF90.Ph+180;
        RFadia.Name = 'sech';
        RFadia.FA = 180;
        RFadia.Ph = PreparationInfo{ii}.RefPh;
        RFadia.OffRes = PreparationInfo{ii}.RefOffRes;
        RFadia.Dura = PreparationInfo{ii}.RefDura;
        RFadia.Npt = 256;
        TM4 = PreparationInfo{ii}.MixTime/4;
        TM2 = PreparationInfo{ii}.MixTime/2;
        [MtxTM4,VecTM4] = func_rotation_decay(paramUse,0,nPool,gamma,TM4);
        [MtxTM2,VecTM2] = func_rotation_decay(paramUse,0,nPool,gamma,TM2);
        Mprep = func_excitation(RF90,paramUse,b1Scale,gamma,Mprep); %90x
        Mprep = MtxTM4*Mprep+VecTM4;
        Mprep = func_excitation(RFadia,paramUse,b1Scale,gamma,Mprep); %adiabatic inversion
        Mprep = MtxTM2*Mprep+VecTM2;
        Mprep = func_excitation(RFadia,paramUse,b1Scale,gamma,Mprep); %adiabatic inversion
        Mprep = MtxTM4*Mprep+VecTM4;
        Mprep = func_excitation(RF90n,paramUse,b1Scale,gamma,Mprep); %90x
    end
end
end