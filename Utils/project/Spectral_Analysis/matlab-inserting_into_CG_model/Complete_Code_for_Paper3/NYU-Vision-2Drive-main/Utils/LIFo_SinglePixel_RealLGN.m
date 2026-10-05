%% LIFonly for Pixels. E is divided into S and C cells
%% Input: 
%        Th idea of this is similar to MF+v, but all L4 input 

%        lgnTF                    LGN temporal frequency
%        L4E, L4I                 Total E/I input kick numbers
%        Assume L4SE/CE/IE are proportional to presynaptic neuron numbers
%        S_EE,S_EI,S_IE,S_II      synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6, L6 Input, All Scalar
%        S_amb,rS_amb,rC_amb,rI_amb      amb Input, All Scalar
%        p_EEFail                 E-to-E Synaptic failure prob
%        lgn_S/COnOff,lgn_I       LGN input rate. Phase OnOff for E (2*1 vec) 
%        S_Elgn S_Ilgn            LGN strength of drives
%        NlgnS,NlgnC,NlgnI        LGN neuron # to SCI cells
%        gL_E,gL_I                Leaky time constants
%        Ve,Vi                    Reversal potentials
%        HyperPara           Cell. Entries are simulation hyperparameters
%                            1th: LIF simulation time
%                            2th an 3th: testmode and test value: threshold
%                            or delay time
%                            
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I

%% Ver 5:       For single pixels
%% Ver 6:       Bg features added:
%% Ver 7:       Using Real LGN data now
%               either lgn_S, lgn_C,lgn_I are bg, then 20 Hz, and N_lgn
%               needed to multiply;
%               or lgn_S and lgn_C are 2*1 cells with LGN spike seqs (S seleted for angle ),
%                  lgn_I is an n*1 cell with multiple LGN spike seqs, and
%                  we propose to randomly select cycles from these seqs
% iteration
% Zhuo-Cheng Xiao 10/23/2023
function [f_EnIOut]...
           = LIFo_SinglePixel_RealLGN(...
...% MF Parameters                     
                     lgnTF, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_S, lgn_C,lgn_I,...%N_Slgn,N_Clgn,N_Ilgn, ...
                     S_Elgn,S_Ilgn,... %7  
...% Below are LIF details
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara) % if varagin non empty, we only export the last state
%% Hyperparameters
N_Hyp = length(HyperPara);

if N_Hyp > 0 % Specify: LIF simulation timef_pre
    LIFSimuT = HyperPara{1};
else
    LIFSimuT = 50*1e3;    % unit in ms    
end

if N_Hyp > 1 % Specify: Starting point of Recording
    if strcmpi(HyperPara{2},'thre')
        RecdThre = HyperPara{3};
        RecdDely = 0;
    elseif strcmpi(HyperPara{2},'delay')
        RecdThre = 0;
        RecdDely = HyperPara{3};
    else
        disp('***Wrong test mode')
        return
    end
else
    RecdThre = 0;    % unit in ms  
    RecdDely = 0;
end
dt = 0.1;

% Note: If rate is per second, need to be divided by 1e3 to rescale to ms
TimeFrac = 1/(LIFSimuT/1e3);%/lgnTF;
    %CycPerSec = floor(LIFSimuT/1e3 * TimeFrac); % should be a whole nubmer, standing for simulation cyc per second.

%% For each cyc, initate a new start.
CellSim = 5; % S4 S6 C1 C2 I4
TCyc = floor(1/TimeFrac); % number of cycs
%VRcdCyc = cell(TCyc,1);
SpRcdCyc = cell(TCyc,1);
RefTimer = zeros(CellSim, 1); % initiate the timer of ref
vt = zeros(CellSim, 1); 
gL = [gL_E*ones(CellSim-1,1); gL_I*ones(1)];
% just for the vector size only
L4Pix_EventsE = PoissonInputForNetwork(CellSim,ones(CellSim,1)/1e3,...
                                       LIFSimuT*TimeFrac,dt,true); % Need to count EE failure here 

% put all lgn_I samples together; and permute. 
if iscell(lgn_S) % drive: use REAL LGN!
    %disp('Drive regime. Use presimulated LGN input')
    TLGNSimu = ceil(max(lgn_S{2})/1e3)*1e3; % should be 24 seconds or 24e3 ms
    Cyctime = 1e3/lgnTF;
    lgn_IAllUse = lgn_I;
    lgn_IAllcyc = cell(length(lgn_IAllUse)*floor(TLGNSimu/1e3)*lgnTF,1);
    %1. first dissect into TF cycs
    cycId = 1;
    for lgn_Iid = 1:length(lgn_IAllUse)
        while cycId <= floor(TLGNSimu/1e3)*lgnTF * lgn_Iid
        lgn_IAllUseNow = lgn_IAllUse{lgn_Iid};
        % start and end time of the cycle
        TStart = (cycId-1)*Cyctime - (lgn_Iid-1)*TLGNSimu; 
        TEnd   =  cycId   *Cyctime - (lgn_Iid-1)*TLGNSimu;
        lgn_IAllUseNow(lgn_IAllUseNow > TEnd | lgn_IAllUseNow < TStart) = [];

        lgn_IAllcyc{cycId} = lgn_IAllUseNow - TStart;
        cycId = cycId + 1;
        end
    end
    %2. then permute the cycles
    lgn_IAllcycPerm = lgn_IAllcyc(randperm(cycId-1));
    for cycId = 1:length(lgn_IAllcycPerm)
        lgn_IAllcycPerm{cycId} = lgn_IAllcycPerm{cycId} + Cyctime * (cycId-1);
    end
    lgn_IAllUseMat = cell2mat(lgn_IAllcycPerm);
end

% Now start L4-LIF simulation
for TInt = 1:TCyc
    VRcd = zeros(size(L4Pix_EventsE)); % this can go away since I don' need to record V-trace
    SpLIF = zeros(CellSim, 1);

    % resample Amb, L6, L4 inputs for every segment;
    % These variables are 1 * Nt vectors
    AmbR = [rS_amb;rS_amb;rC_amb;rC_amb;rI_amb]; % rQ_amb are per ms
    L6R  = [rL6SU; rL6SU; rL6CU; rL6CU; rL6IU];  % rL6 are per ms
    L4ER = [L4SEU*(1-p_EEFail); 
            L4SEU*(1-p_EEFail); 
            L4CEU*(1-p_EEFail); 
            L4CEU*(1-p_EEFail); 
            L4IEU]/1e3; % Need to count EE failure here 
    L4IR = [L4SIU; 
            L4SIU; 
            L4CIU; 
            L4CIU; 
            L4IIU]/1e3; % L4E/I are per second
    AmbPix_Events = PoissonInputForNetwork(CellSim,AmbR,LIFSimuT*TimeFrac,dt,true); 
    L6Pix_Events  = PoissonInputForNetwork(CellSim,L6R, LIFSimuT*TimeFrac,dt,true); 
    L4Pix_EventsE = PoissonInputForNetwork(CellSim,L4ER,LIFSimuT*TimeFrac,dt,true); 
    L4Pix_EventsI = PoissonInputForNetwork(CellSim,L4IR,LIFSimuT*TimeFrac,dt,true); 

    % LGN input pissons are handeled specifically:
    if iscell(lgn_S) % drive: use REAL LGN!
        %TLGNSimu = ceil(max(lgn_S{2})/1e3)*1e3; % should be 24 seconds
        % get the corresponding segment used from LGN simulations
        lgnEvents = cell(CellSim,1);

        TStart = (TInt-1)*LIFSimuT*TimeFrac; TEnd = TInt*LIFSimuT*TimeFrac;
        TStartE = mod(TStart,TLGNSimu);       TEndE = mod(TEnd,TLGNSimu);
        if TEndE == 0
            TEndE = TLGNSimu;
        end
        % need to re-align
        lgnEvents{1} = lgn_S{1}(lgn_S{1}>TStartE & lgn_S{1}<TEndE) - TStartE; % lgnS4
        lgnEvents{2} = lgn_S{2}(lgn_S{2}>TStartE & lgn_S{2}<TEndE) - TStartE; % lgnS6
        lgnEvents{3} = lgn_C{1}(lgn_C{1}>TStartE & lgn_C{1}<TEndE) - TStartE; % lgnC1
        lgnEvents{4} = lgn_C{2}(lgn_C{2}>TStartE & lgn_C{2}<TEndE) - TStartE; % lgnC2
        lgnEvents{5} = lgn_IAllUseMat(lgn_IAllUseMat>TStart & lgn_IAllUseMat<TEnd) - TStart; % lgnI4
        
        lgn_Events = zeros(size(AmbPix_Events));
        for LIFId = 1:CellSim
            lgn_Events(LIFId,ceil(lgnEvents{LIFId}/dt)) = 1;
        end
    elseif isnumeric(lgn_S) && length(lgn_S) == 1 % bg
        %disp('Background regime. Use Phase consistent LGN input')
        lgnBGR = [4;6;1;2;4].*[lgn_S;lgn_S;lgn_C;lgn_C;lgn_I]; % first vec N_lgn used:S4 S6 C1 C2 I4
        lgn_Events = PoissonInputForNetwork(CellSim,lgnBGR,LIFSimuT*TimeFrac,dt,true); % lgn_Q should be per ms
    else
        disp('***Illigal LGN input. Returning...')
    end

    lgnPix_Events  = lgn_Events;% should pack up all LGN and L4
    % New Gs
    winAMPA = 0:dt:tau_ampa_D*3;
    winNMDA = 0:dt:tau_nmda_D*3;
    winGABA = 0:dt:tau_gaba_D*3;
    KerAMPA = 1/(tau_ampa_D-tau_ampa_R) * (exp(-winAMPA/tau_ampa_D) - exp(-winAMPA/tau_ampa_R));
    KerNMDA = 1/(tau_nmda_D-tau_nmda_R) * (exp(-winNMDA/tau_nmda_D) - exp(-winNMDA/tau_nmda_R));
    KerGABA = 1/(tau_gaba_D-tau_gaba_R) * (exp(-winGABA/tau_gaba_D) - exp(-winGABA/tau_gaba_R));
    KerAMPA = KerAMPA / (sum(KerAMPA)*dt);
    KerNMDA = KerNMDA / (sum(KerNMDA)*dt);
    KerGABA = KerGABA / (sum(KerGABA)*dt);
    % Incorporate L4 and other input
    ampaInp = single([S_Elgn * lgnPix_Events(1:CellSim-1,:);            S_Ilgn * lgnPix_Events(end,:)]...
                    +  S_amb  * AmbPix_Events ... % S_amb identical for E and I
                    + [S_EL6 * L6Pix_Events(1:CellSim-1,:) * rhoE_ampa; S_IL6 * L6Pix_Events(end,:) * rhoI_ampa]...
                    + [S_EE * L4Pix_EventsE(1:CellSim-1,:) * rhoE_ampa; S_IE * L4Pix_EventsE(end,:) * rhoI_ampa]);
    nmdaInp = single([S_EL6  * L6Pix_Events(1:CellSim-1,:) * rhoE_nmda; S_IL6 * L6Pix_Events(end,:) * rhoI_nmda]...
                    + [S_EE * L4Pix_EventsE(1:CellSim-1,:) * rhoE_nmda; S_IE * L4Pix_EventsE(end,:) * rhoI_nmda]);
    gabaInp = single([S_EI * L4Pix_EventsI(1:CellSim-1,:);              S_II * L4Pix_EventsI(end,:)] );
    %% Precompute All Conductances: Speed up by ifft/fft
    n = size(ampaInp,2);
    GAMPA = ifft(fft(ampaInp',n) .* repmat(fft(KerAMPA',n),1,CellSim))';
    GNMDA = ifft(fft(nmdaInp',n) .* repmat(fft(KerNMDA',n),1,CellSim))';
    GGABA = ifft(fft(gabaInp',n) .* repmat(fft(KerGABA',n),1,CellSim))';

    GE = GAMPA + GNMDA;  GI = GGABA;
    
    for tInd = 1:size(L4Pix_EventsE,2)
        % Firstly, refrectory neurons get out if timer reachs t_ref
        RefTimer(isnan(vt)) = RefTimer(isnan(vt)) + dt;
        vt(RefTimer>=tau_ref) = 0;
        RefTimer(RefTimer>=tau_ref) = 0;

        ov  = vt + dt*(-gL.*vt - GE(:,tInd).*(vt-Ve) - GI(:,tInd).*(vt-Vi));
        SpNow = single(ov>=1);
        ov(ov>=1) = nan;
        VRcd(:,tInd) = ov;
        % Compute meanV
            SpLIF = SpLIF + SpNow;
        % Transfer variables
        vt = ov;
    end
    %VRcdCyc{TInt} = VRcd;
    SpRcdCyc{TInt} = SpLIF;
end

%% I am not concerning the Vs now, but maybe need that in future
%VRcdMat = cell2mat(VRcdCyc')';
SpRcdMat = cell2mat(SpRcdCyc')';
f_EnIOut = SpRcdMat'*(1e3/(LIFSimuT*TimeFrac)); 
                    % 1e3/(LIFSimuT*TimeFrac): this is how many cycles per
                    % second
%FrRcdAve(LIFtestInd,:) = mean(SpRcdMat(10:end,:))*lgnTF;


end





