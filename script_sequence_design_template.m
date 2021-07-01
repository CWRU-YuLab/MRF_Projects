%% Input experimental setup
clear all
b1Scale = 1; %b1 scale, actual flip angle = designed flip angle * b1Scale
gamma = 17.235; %gyromagnetic ratio for P31, unit: MHz/T
% gamma = 42.58; %gyromagnetic ratio for proton, unit: MHz/T
b0 = 9.4; %main field strength, unit: T

%% -------- Step 1: Generate the sequence parameters -----------
%  This part has hard-coded sequence parameters
%  Save the script for bookkeeping purposes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------- General sequence design -----------------
% Put down general description of the sequence here for book keeping
%
% Example: 3-segment P31-MRSF sequence
%
% Segment 1: global inversion, alternating PCr and ATP acquisition
% 8 blocks of PCr and 8 blocks of ATP, use the same flip angle pattern
%
% Segment 2: ATP saturation, interleaved with PCr acquisition
% 8 blocks only, use the same flip angle pattern
%
% Segment 3: PCr/Pi saturation, interleaved with ATP acquisition
% 8 blocks only, use the same flip angle pattern

NR = 320; % Total number of TRs in the sequence
N_segments = 3; % Total number of segments
N_blocks = zeros(1,N_segments);
N_blocks(1) = 16; % Number of acquisition blocks in Segment 1
N_blocks(2) = 8; % Number of acquisition blocks in Segment 2
N_blocks(3) = 8; % Number of acquisition blocks in Segment 3

load('FA_31P_MRSF'); % Load flip angle files

%-------------------- Set the RF pulse in each TR -------------------------
% GaussTrunc20 for all excitation pulses, 
% Pulse duration is 4-ms, simulated in 41 points
% FA pattern predesigned and loaded from a file
% Alternating phase between 0 and 180

excName = repmat({'GaussTrunc20'},1,NR);
excFA = reshape(FA,1,NR); % Designed flip angle, unit: degree
excPh = zeros(1,NR); % RF phase relative to the x-axis of the rotating frame, unit: degree
excPh(2:2:end) = 180;
excOffRes = zeros(1,NR); % RF off-resonance frequency to the rotating frame, unit: Hz
for ii = 1:N_segments
    switch ii
        case 1 % segment 1, 16 blocks alternating between PCr and ATP, 10 TRs each
            for jj = 1:2:16
                excOffRes(1,(10*jj+1):10*(jj+1)) = gamma*b0*(2.4); %ATP excitation
            end
        case 2 % segment 2, 8 blocks, PCr only
        case 3 % segment 3, 8 blocks, ATP only, from TR# 241 to NR
            excOffRes(1,241:NR) = gamma*b0*(2.4);
    end
end
excDura = 0.004*ones(1,NR); %RF duration, unit: s
excNpt = 41*ones(1,NR); %Number of data points to simulate the RF shape

%-------- Set the timging of data acquisition in each TR -------------
TR = 0.01*ones(1,NR); % Constant repetition time of 0.01 s
adcDelay = 0.000001*ones(1,NR); % Delay between end of RF pulse and start of ADC, unit: s
adcDura = 0.007725*ones(1,NR); % Readout duration, unit: s
adcNpt =  70*ones(1,NR); % Number of data points in one readout
dePhase = 0*ones(1,NR); % Dephase related to net gradient within one TR

%-------- Set preparation modules ------------
% ########### Pre-defined praparation modules, DO NOT DELETE #############
% Set the parameters before adding the module to the sequence

% Inversion, 1-ms duration, 10-us delay
inv = struct('Type','Inversion','InvName','bp','InvFA',180,...
    'InvPh',0,'InvOffRes',0,'InvDura',0.001,'InvNpt',1,'InvDelay',0.00001);

% Adiabatic inversion, 10-ms duration, 10-us delay
invSech = struct('Type','InversionSech','InvFA',180,...
    'InvPh',0,'InvOffRes',0,'InvDura',0.01,'InvDelay',0.00001);

% Spectral selective ideal inversion
invSpec = struct('Type','InversionSpecSel','PoolID',[1 2 3]);

% Continuous-wave saturation, 485-ms duration
% SatB1 unit:uT, current value = rot 5 degree/0.5ms
satBp = struct('Type','Bp-Saturation','SatB1',5/180*pi/(2*pi*gamma*0.0005),...
    'SatPh',0,'SatOffRes',0,'SatDura',0.485); 

% Spectral selective ideal saturation, PCr/Pi saturation, 485-ms duration
satSpec = struct('Type','SpecSelSaturation','PoolID',[1 3],'SatDura',0.485);

% Crusher, 5-ms duration
crusher = struct('Type','Crusher','CrushDura',0.005);

% Delay, 0.5-s duration
delay = struct('Type','Delay','DelayDura',0.5);

% T2 prep with MLEV4 refocusing
t2_MLEV4 = struct('Type','T2-MLEV4','ExcPh',0,'ExcOffRes',0,...
    'ExcDura',0.0003,'MixTime',0.02);

% T2 prep with adiabatic refocusing
t2_Adia = struct('Type','T2-Adiabatic','ExcPh',0,'ExcOffRes',0,'ExcDura',0.0003,...
    'RefPh',0,'RefOffRes',0,'RefDura',0.008,'MixTime',0.02);

% ################ End of Pre-Defined Praparation Modules ################

% Arrange prep modules in MRF sequence 
prepMod = cell(1,NR+1); % Setup NR+1 array to accommodate a delay at the end

% Global inversion at the beginning, followed by crusher
invSpec.PoolID = [1 2 3]; 
prepMod{1} = cell(1,1);
prepMod{1}{1} = invSpec;
prepMod{1}{2} = crusher;

% Control saturation in Segment 1
satBp.SatOffRes = -gamma*b0*2.4; % control saturation
satSpec.SatDura = 0.485; % 485-ms saturation
for ii = 11:10:161
    prepMod{ii} = cell(1,1);
    prepMod{ii}{1} = satBp;
    prepMod{ii}{2} = crusher;
end

% ATP saturation in Segment 2
satBp.SatOffRes = gamma*b0*2.4; % gATP saturation
satSpec.SatDura = 0.485; % 485-ms saturation
satSpec.PoolID = [2]; %ATP saturation
for ii = 171:10:241
    prepMod{ii} = cell(1,1);
    prepMod{ii}{1} = satSpec;
    prepMod{ii}{2} = crusher;
end

% PCr&Pi saturation in Segmentation 3 
satSpec.SatDura = 0.485; % 485-ms saturation
satSpec.PoolID = [1 3]; %PCr&Pi saturation
for ii = 251:10:320
    prepMod{ii} = cell(1,1);
    prepMod{ii}{1} = satSpec;
    prepMod{ii}{2} = crusher;
end


%% --------- Step 2: Assemble sequence into the data structure "Sequence" ------------
clear Sequence

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
Sequence.NR = NR;
gtext = input('Enter a short note:','s');
Sequence.Note = gtext;
uisave({'Sequence'})
