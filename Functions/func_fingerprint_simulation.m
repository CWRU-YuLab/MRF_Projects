function Fingerprint = func_fingerprint_simulation(Sequence, Variable, b1Scale, gamma, b0, dummy)
% This function simulates a single fingerprint using sequence and tissue
% parameters as inputs
%
% Inputs:
%   Sequence: contains all parameters of an MRF sequences in each TR, including
%             preparation modules preceding each TR
%   Variable: contains all tissue parameters needed for fingerprint simulation
%             (input variable is repackaged into tissue_parameter using
%             func_parameter_generation
% 
% Outputs:
%   entry: fingerprint entries, dim1-Mxyz1Mxyz2... (3*#compartment)
%                      dim2-#pts in FID, dim3-#TR

%findout all parameters
nPool = length(Variable.pool);
NR = length(Sequence.RF);
tissue_parameters = func_parameter_generation(Variable,1,1);
n_tissue_parameters = length(tissue_parameters); % this is necessary!! -Kihwan Kim- 

%loop through all dictionary entries
adcNpt = max([Sequence.FP.adcNpt]);
offResApp = tissue_parameters.OffRes;%seqDMpar(ii,(nPool*3+1):(nPool*4));
offResAll = repmat(reshape(linspace(Variable.vOffResMax(1),Variable.vOffResMax(2),Variable.nSpins),[],1),1,nPool);
offResAll = offResAll+repmat(offResApp,Variable.nSpins,1);
seqFid = zeros(nPool*3,adcNpt,NR,Variable.nSpins);

for kk = 1:n_tissue_parameters;
    paramUse = tissue_parameters(kk); % loop through the number of tissue parameter entries and assign 
    M0 = zeros(nPool*3,1);
    M0(3:3:end) = paramUse.M0;
    for jj = 1:Variable.nSpins
        paramUse.OffRes = offResAll(jj,:);
        Mprev = M0;
        for dd = dummy:-1:0
            for kk = 1:NR
                %preparation
                if ~isempty(Sequence.Prep(kk))
                    Mprep = func_preparation(Sequence.Prep(kk).Module,paramUse,b1Scale,gamma,Mprev);
                else
                    Mprep = Mprev;
                end
                
                %excitation
                RF = Sequence.RF(kk);
                [Mrf Jacob] = func_excitation(RF,paramUse,b1Scale,gamma,Mprep);
                
                % Free Precession and acquisition
                % timing of pre adc delay
                FP = Sequence.FP(kk);
                if dd~=0 %not dummy scan
                    FP.adcDura = 0;%not used
                    FP.adcNpt = 1;
                end
                %timging of dephase duration
                FP.dePhaseDura = FP.TR-RF.Dura-FP.adcDelay;
                %off-res frequency of spin jj caused by gradient
                if(FP.dePhase~=0)
                    dePhaseMax = FP.dePhase/2/pi/FP.dePhaseDura;
                    FP.curDePhase = (jj-1)/(Variable.nSpins-1)*dePhaseMax;
                    FP.curDePhase = FP.curDePhase - dePhaseMax/2;
                else
                    FP.curDePhase = 0;
                end
                %calculate magnetization after RF to end of TR
                [Madc,Mtr] = func_freeprecession(FP,paramUse,Mrf);
                
                if dd==0
                    fidTemp = Madc.*Variable.wVOffRes(jj);
                    %readout phase correction
                    for nn = 1:nPool
                        temp = (fidTemp(1+(nn-1)*3,:)+1i*fidTemp(2+(nn-1)*3,:))*exp(1i*RF.Ph/180*pi); %ADC phase correction
                        fidTemp(1+(nn-1)*3,:) = real(temp);
                        fidTemp(2+(nn-1)*3,:) = imag(temp);
                    end
                    seqFid(:,:,kk,jj) = fidTemp;
                end
                Mprev = Mtr;
                
            end
            if ~isempty(Sequence.Prep(NR+1))
                Mprev = func_preparation(Sequence.Prep(NR+1).Module,paramUse,b1Scale,gamma,Mprev);
            end
        end
    end
end
Fingerprint = paramUse;
Fingerprint.Entry = squeeze(sum(seqFid,4));