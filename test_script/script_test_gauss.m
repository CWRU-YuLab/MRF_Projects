%% Input experimental setup
b1Scale = 1; %b1 scale, actual flip angle = designed flip angle * b1Scale
gamma = 17.23; %gyromagnetic ratio, unit: MHz/T
b0 = 9.4; %main field strength, unit: T

%% Imaging sequence parameters
clear Sequence
%-------------------- Set the RF pulse in each TR -------------------------
NR = 1; %Number of TR in the sequence
excName = repmat({'GaussTrunc20'},1,NR);
excFA = 90*ones(1,NR); %Designed flip angle, unit: degree
excPh = zeros(1,NR); %RF phase relative to the x-axis of the rotating frame, unit: degree
excOffRes = 2.4*gamma*b0*ones(1,NR); %RF off-resonance frequency to the rotating frame, unit: Hz
excDura = 0.004*ones(1,NR); %RF duration, unit: s
excNpt = 501*ones(1,NR); %Number of data points to simulate the RF shape

%-------- Set the timging of data acquisition in each TR -------------
adcDelay = 0.000001*ones(1,NR); %Delay between end of RF pulse and start of ADC, unit: s
adcDura = 0.007725*ones(1,NR); %Readout duration, unit: s
adcNpt =  70*ones(1,NR); %Number of data points in one readout
dePhase = 0*ones(1,NR); %Dephase related to net gradient within one TR
TR = 0.01*ones(1,NR);%repetition time of the imaging sequence

%-------- Set prep modules before each acquisition module ------------
%structures predefined for different prep modules
%inversion
inv = struct('Type','Inversion','InvName','bp','InvFA',180,...
    'InvPh',0,'InvOffRes',0,'InvDura',0.001,'InvNpt',1,'InvDelay',0.00001);
%adiabatic inversion
invSech = struct('Type','InversionSech','InvFA',180,...
    'InvPh',0,'InvOffRes',0,'InvDura',0.01,'InvDelay',0.00001);
%continuous-wave saturation
satBp = struct('Type','Bp-Saturation','SatB1',5/180*pi/(2*pi*gamma*0.0005),...
    'SatPh',0,'SatOffRes',0,'SatDura',0.485); %sat B1 unit:uT, current value = rot 5 degree/0.5ms
%crusher
crusher = struct('Type','Crusher','CrushDura',0.005);
%delay
delay = struct('Type','Delay','DelayDura',0.5);
%T2 prep with MLEV4 refocusing
t2_MLEV4 = struct('Type','T2-MLEV4','ExcPh',0,'ExcOffRes',0,...
    'ExcDura',0.0003,'MixTime',0.02);
%T2 prep with adiabatic refocusing
t2_Adia = struct('Type','T2-Adiabatic','ExcPh',0,'ExcOffRes',0,'ExcDura',0.0003,...
    'RefPh',0,'RefOffRes',0,'RefDura',0.008,'MixTime',0.02);
%arrange prep modules in MRF sequence 
prepMod = cell(1,NR+1); %+1 because probably a delay can be inserted after NR

%%
%--------- Pass sequence design into data structure "Sequence" ------------
Sequence = struct('Prep',[],'RF',[],'FP',[]);
%prepMod is 1x(NR+1) cell matrix
Sequence.Prep = struct('Module',prepMod);
%excName is 1xNR cell matrix, excFA/Ph/OffRes/Dura/Npt are 1xNR vectors
Sequence.RF = struct('Name',excName, 'FA',num2cell(excFA,1),...
    'Ph',num2cell(excPh,1),'OffRes',num2cell(excOffRes,1),...
    'Dura',num2cell(excDura,1),'Npt',num2cell(excNpt,1));
%adcDelay/Dura/Npt, dePhase, TR are 1xNR vectors
Sequence.FP = struct('adcDelay',num2cell(adcDelay,1),'adcDura',num2cell(adcDura,1),...
    'adcNpt',num2cell(adcNpt,1),'dePhase',num2cell(dePhase,1),...
    'TR',num2cell(TR,1));


%% Generate model variables for dictionary simulation
clear Variable
%input the range and step size of variables here
%------------- Set T1, T2, M0, off-res for each compartment ---------------
Variable.pool = struct('T1',[],'T2',[],'M0',[],'OffRes',[],...
    'rOffResPool',[],'rOffRes',[],'OutFluxRate',[]);
nPool = 2;
%for example, pool 1 is assigned for PCr
Variable.pool(1).T1 = 3.3; %T1,unit: s
Variable.pool(1).T2 = 0.12; %T2,unit: s
Variable.pool(1).M0 = 3; %equilibrium magnetization, a.u.
Variable.pool(1).OffRes = 0; %apparent off resonance in one voxel, relative to rotating frame, unit: Hz
Variable.pool(1).OutFluxRate = cell(1,nPool);
Variable.pool(1).OutFluxRate{2} = 0; %rate of MT from PCr to ATP
%for example, pool 2 is assigned for ATP 
Variable.pool(2).T1 = 0.8;
Variable.pool(2).T2 = 0.016; 
Variable.pool(2).M0 = 1;
Variable.pool(2).OffRes = []; %the off res will be a fixed value relative another pool
Variable.pool(2).rOffResPool = 1; %offset will be relative to pool1
Variable.pool(2).rOffRes = gamma*b0*2.4; %offset value (Hz)
Variable.pool(2).OutFluxRate = cell(1,nPool);
Variable.pool(2).OutFluxRate{1} = 'Equil'; %flux of MT from ATP to PCr will be in equilibrium with PCr to ATP 

%-------- Set off-resonance within voxel shared by all compartments -------
Variable.nSpins = 1; %number of spins in one voxel 
Variable.vOffResMax = [0; 0];%maximal off res relative to the apparent off res in one voxel, unit: Hz
Variable.wVOffRes = 1; %the probability of each off res value 
Variable.wVOffRes = Variable.wVOffRes./sum(Variable.wVOffRes);

%%
nSubDic = 1;
nDic = 1;
dummy = 0;

tic
DictionaryFID = func_dictionary_generation(Sequence, Variable, b1Scale, gamma, b0, dummy, nSubDic, nDic);
toc

%FID should be sum of Mxy in PCr and ATP 
fidAll = cat(4,DictionaryFID.Entry); %extract all dictionary entries
entryFIDsum = sum(fidAll(1:3:end,:,:,:),1)+1i*sum(fidAll(2:3:end,:,:,:),1); %Mx+1i*My
% do FFT along time domain to get spectrum
entrySpec = fftshift(fft(ifftshift(entryFIDsum,2),[],2),2);
freqMax = 1/(adcDura(1)/adcNpt(1))/2;
freq = linspace(freqMax,-freqMax,adcNpt(1))/gamma/b0;

%assign fingerprint to dictionary structure
Dictionary = struct('T1',{DictionaryFID.T1},'T2',{DictionaryFID.T2},...
    'M0',{DictionaryFID.M0},'OffRes',{DictionaryFID.OffRes},...
    'MtRate',{DictionaryFID.MtRate},...
    'Entry',squeeze(mat2cell(entrySpec,1,adcNpt(1),NR,ones(1,length(DictionaryFID))))');
%%
figure(1)
hold on
plot(repmat(freq,[length(Dictionary) 1])',abs(cat(1,Dictionary.Entry))','LineWidth',1.5)
xlabel('Chemical Shift (ppm)')
ylabel('Magnitude (a.u.)')
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'box','on')
xlim([-15 15])