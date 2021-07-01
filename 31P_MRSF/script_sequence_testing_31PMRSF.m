%% Input experimental setup
clear all
b1Scale = 1; %b1 scale, actual flip angle = designed flip angle * b1Scale
gamma = 17.235; %gyromagnetic ratio, unit: MHz/T
b0 = 9.4; %main field strength, unit: T

uiload;
%% Generate model variables for dictionary simulation (2 pool)
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
Variable.pool(1).OutFluxRate{2} = [0 0.5]; %rate of MT from PCr to ATP

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
Variable.nSpins = 51; %number of spins in one voxel 
Variable.vOffResMax = [-75; +75];%maximal off res relative to the apparent off res in one voxel, unit: Hz
lwfwhm = 15; %prepare parameter for Lorentzian linewidth, 1 value at a time, unit: Hz
Variable.wVOffRes = func_lineshape('Lorentzian',Variable.nSpins,lwfwhm); %the probability of each off res value 
Variable.wVOffRes = Variable.wVOffRes./sum(Variable.wVOffRes);
%% Generate model variables for dictionary simulation (3 pool)
clear Variable
%input the range and step size of variables here
%------------- Set T1, T2, M0, off-res for each compartment ---------------
Variable.pool = struct('T1',[],'T2',[],'M0',[],'OffRes',[],...
    'rOffResPool',[],'rOffRes',[],'OutFluxRate',[]);
nPool = 3;
%for example, pool 1 is assigned for PCr
Variable.pool(1).T1 = 4.0; %T1,unit: s
Variable.pool(1).T2 = 0.22; %T2,unit: s
Variable.pool(1).M0 = 3; %equilibrium magnetization, a.u.
Variable.pool(1).OffRes = 0; %apparent off resonance in one voxel, relative to rotating frame, unit: Hz
Variable.pool(1).OutFluxRate = cell(1,nPool);
Variable.pool(1).OutFluxRate{2} =[0 0.5]; %rate of MT from PCr to ATP

%for example, pool 2 is assigned for ATP 
Variable.pool(2).T1 = 3.3;
Variable.pool(2).T2 = 0.029; 
Variable.pool(2).M0 = 1;
Variable.pool(2).OffRes = []; %the off res will be a fixed value relative another pool
Variable.pool(2).rOffResPool = 1; %offset will be relative to pool1
Variable.pool(2).rOffRes = gamma*b0*2.4; %offset value (Hz)
Variable.pool(2).OutFluxRate = cell(1,nPool);
Variable.pool(2).OutFluxRate{1} = 'Equil'; %flux of MT from ATP to PCr will be in equilibrium with PCr to ATP 
Variable.pool(2).OutFluxRate{3} = 'Equil';

%for example, pool 2 is assigned to Pi
Variable.pool(3).T1 = 6.3;
Variable.pool(3).T2 = 0.109; 
Variable.pool(3).M0 = 1;
Variable.pool(3).OffRes = []; %the off res will be a fixed value relative another pool
Variable.pool(3).rOffResPool = 1; %offset will be relative to pool1
Variable.pool(3).rOffRes = gamma*b0*(-5.0); %offset value (Hz)
Variable.pool(3).OutFluxRate = cell(1,nPool);
Variable.pool(3).OutFluxRate{2} = 0.1; %flux of MT from ATP to PCr will be in equilibrium 

%-------- Set off-resonance within voxel shared by all compartments -------
Variable.nSpins = 51; %number of spins in one voxel 
Variable.vOffResMax = [-75; +75];%maximal off res relative to the apparent off res in one voxel, unit: Hz
lwfwhm = 15; %prepare parameter for Lorentzian linewidth, 1 value at a time, unit: Hz
Variable.wVOffRes = func_lineshape('Lorentzian',Variable.nSpins,lwfwhm); %the probability of each off res value 
Variable.wVOffRes = Variable.wVOffRes./sum(Variable.wVOffRes);

%% Perform Bloch simulation
dummy = 1;
%seqDMpar: Number of dictionary entry, Number of Variables
%seqDMdicFID: 3*nPool(Mx1,My1,Mz1,...), adcNpt, NR, Number of Variables  
tic
fingerprint = func_fingerprint_simulation(Sequence, Variable, b1Scale, gamma, b0, dummy);
toc
%%
%FID should be sum of Mxy in PCr and ATP 
fidAll = cat(4,fingerprint.Entry); %extract all dictionary entries
entryFIDsum = squeeze(sum(fidAll(1:3:end,:,:,:),1)+1i*sum(fidAll(2:3:end,:,:,:),1)); %Mx+1i*My
% do FFT along time domain to get spectrum
entrySpec = squeeze(fftshift(fft(ifftshift(entryFIDsum,1),[],1),1));
% extract signal corresponding to a metabolite
% PCr: 36
% ATP: 33
% Pi:  42
pcrSignal = reshape(multiaa(squeeze(entrySpec(36,:,:))),Sequence.NR,[]);
atpSignal = reshape(multiaa(squeeze(entrySpec(42,:,:))),Sequence.NR,[]);
fingerprint.PCr = pcrSignal;
fingerprint.ATP = atpSignal;
gtext = input('Enter a short note:','s');
fingerprint.Note = gtext;
uisave({'fingerprint'})
%%
figure
subplot(2,1,1)
plot(real(pcrSignal),'LineWidth',1.5); hold on
plot(imag(pcrSignal),'LineWidth',1.5);
xlim([0 320])
set(gca,'XTick',0:80:320)
ylabel('PCr(a.u.)')
xlabel('#TR')
gtext=sprintf('k^f_C_K=%1.2f, k^f_A_T_P=%1.2f',fingerprint.MtRate(1,2),fingerprint.MtRate(3,2));
title(gtext)
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'box','on')
subplot(2,1,2)
plot(real(atpSignal),'LineWidth',1.5); hold on
plot(imag(atpSignal),'LineWidth',1.5); 
xlim([0 320])
%ylim([-80 120])
set(gca,'XTick',0:80:320)
ylabel('Pi(a.u.)')
% legend('kf=0','kf=0.5','kf=0','kf=0.5')
xlabel('#TR')
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'box','on')