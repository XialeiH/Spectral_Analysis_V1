function simulate_paper3_snn_pixels()
% Run one Paper3 SNN condition and save compact pixelwise firing-rate data.

projectRoot = getenv('PAPER3_PROJECT_ROOT');
if isempty(projectRoot)
    projectRoot = '/scratch/xh2906/NYU-Vision-2Drive-main';
end
outputDir = getenv('SNN_OUTPUT_DIR');
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'Data', 'Paper2_NetworkTuning', ...
        'Fig1V4', 'Paper3NWSimulationData_h96_validation');
end

angleDeg = envNumber('SNN_ANGLE', 0);
replicate = round(envNumber('SNN_REPLICATE', 1));
durationMs = envNumber('SNN_DURATION_MS', 10000);
analysisMs = envNumber('SNN_ANALYSIS_MS', min(9000, durationMs));
if analysisMs <= 0 || analysisMs > durationMs
    error('SNN_ANALYSIS_MS must be in (0, SNN_DURATION_MS].');
end

seed = 910000 + round(100 * angleDeg) * 10 + replicate;
rng(seed, 'twister');
addpath(projectRoot, '-begin');
addpath(fullfile(projectRoot, 'Utils'), '-begin');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

dataRoot = fullfile(projectRoot, 'Data', 'Paper2_NetworkTuning');
workspaceFile = fullfile(dataRoot, 'Fig1V4', 'AllMFPixPara_Paper2TuneFig1V4D2.mat');
connectionFile = fullfile(dataRoot, 'Fig1V4', 'Paper3_LIF_ConnMats_V1.mat');
assert(isfile(workspaceFile), 'Missing Paper3 workspace: %s', workspaceFile);
assert(isfile(connectionFile), 'Missing Paper3 connectivity: %s', connectionFile);

W = load(workspaceFile);
C = load(connectionFile, 'C_EE_Fix_Bd', 'C_EI_Fix_Bd', ...
    'C_IE_Fix_Bd', 'C_II_Fix_Bd');
C_EE = C.C_EE_Fix_Bd;
C_EI = C.C_EI_Fix_Bd;
C_IE = C.C_IE_Fix_Bd;
C_II = C.C_II_Fix_Bd;
EcplxInd = logical(W.EcplxInd(:));

N_E = W.N_E; N_I = W.N_I; N_S = W.N_S; N_C = W.N_C;
N_HC = W.N_HC; NPixX = W.NPixX; NPixY = W.NPixY;
dt = 0.1;
TPar = 10;
stepsTotal = round(durationMs / dt);
if mod(stepsTotal, TPar) ~= 0
    error('SNN_DURATION_MS/dt must be divisible by %d.', TPar);
end
stepsPerSection = stepsTotal / TPar;

% Orientation-dependent LGN and L6 drives from
% Paper3_TestL6ShapeDriveHPC_Final.m.
ODNum = 4;
gratingHC = angleDeg + [0, 45, 90, 135];
scaledLGN = (cosd(abs(mod(gratingHC, 180)) * 2) + 1) / 2;
phaseFRS = [45 + scaledLGN * 45, 45 - scaledLGN * 45] * W.N_Slgn / 1e3;
phaseFRC = 45 * ones(1, 8) * W.N_Clgn / 1e3;
OD_E = [W.OD_S, W.OD_C];

L6up = 60; L6low = 6;
L6OrdFHz = L6low + (L6up - L6low) * L6CurveFinalize(scaledLGN);
L6SDrive = L6OrdFHz * W.NS_L6 / 1e3;
L6CDrive = L6OrdFHz * W.NC_L6 / 1e3;
L6IDrive = L6OrdFHz * W.NI_L6 / 1e3;
L6smear = W.n_S_HC * 0.33;
L6trunc = 0.33 / (L6smear / W.n_S_HC);
rS_L6 = SpatialGaussianFilt_my(W.OD_SMap, N_HC, W.n_S_HC, ...
    L6smear, L6trunc, false, L6SDrive);
rC_L6 = SpatialGaussianFilt_my(W.OD_CMap, N_HC, W.n_C_HC, ...
    L6smear, L6trunc, false, L6CDrive);
rI_L6 = SpatialGaussianFilt_my(W.OD_IMap, N_HC, W.n_I_HC, ...
    L6smear, L6trunc, false, L6IDrive);
rE_L6 = [rS_L6; rC_L6];

lambdaEPre = 45 * W.N_Slgn / 1e3;
lambdaIPre = 45 * W.N_Ilgn / 1e3;
inputPoolMs = durationMs * 0.05;
lgnEEvents = PoissonInputForNetwork(N_E, lambdaEPre, inputPoolMs, dt);
lgnIEvents = PoissonInputForNetwork(N_I, lambdaIPre, inputPoolMs, dt);
ambEEvents = PoissonInputForNetwork(N_E, W.rE_amb, inputPoolMs, dt);
ambIEvents = PoissonInputForNetwork(N_I, W.rI_amb, inputPoolMs, dt);
l6EEvents = PoissonInputForNetwork(N_E, rE_L6, inputPoolMs, dt);
l6IEvents = PoissonInputForNetwork(N_I, rI_L6, inputPoolMs, dt);
maxInputBin = floor(inputPoolMs / dt);

% Use the compatible saved Paper2 state, matching the source initialization.
state = W.EndState;
RefTimeE = state{1}; VE = state{2}; SpE = state{3};
GE_ampa_R = state{4}; GE_nmda_R = state{5}; GE_gaba_R = state{6};
GE_ampa_D = state{7}; GE_nmda_D = state{8}; GE_gaba_D = state{9};
RefTimeI = state{10}; VI = state{11}; SpI = state{12};
GI_ampa_R = state{13}; GI_nmda_R = state{14}; GI_gaba_R = state{15};
GI_ampa_D = state{16}; GI_nmda_D = state{17}; GI_gaba_D = state{18};
EEDlyRcd = sparse(N_E, 1);
IEDlyRcd = ones(N_I, 1);
EIDlyRcd = ones(N_E, 1);

LGNFreq = 4;
tMod = 1e3 / LGNFreq;
maxOnPhase = tMod / 2;
PhaseE = single(tMod * rand(N_E, 1));
frameMs = 100;
frameCount = round(frameMs / dt);
tt = 0:dt:(frameMs - dt);
spikeCountE = zeros(N_E, 1);
spikeCountI = zeros(N_I, 1);
analysisStartMs = durationMs - analysisMs;
timerAll = tic;

for section = 1:TPar
    timerSection = tic;
    for timeN = 1:stepsPerSection
        PhaseE = mod(PhaseE + dt, tMod);
        frameIndex = mod(timeN, frameCount);
        if frameIndex == 0; frameIndex = frameCount; end
        if frameIndex == 1
            randomBins = randi([1, maxInputBin + 1], 6, frameCount);
            phaseIndex = floor(mod(PhaseE + tt, tMod) / maxOnPhase) * ODNum + OD_E';
            adjustedLGNE = phaseFRS(phaseIndex);
            adjustedLGNE(EcplxInd,:) = phaseFRC(phaseIndex(EcplxInd,:));
            EampaInp = single(full( ...
                W.S_Elgn * (lgnEEvents(:, randomBins(1,:)) .* adjustedLGNE / lambdaEPre) + ...
                W.S_amb * ambEEvents(:, randomBins(3,:)) + ...
                W.S_EL6 * l6EEvents(:, randomBins(5,:)) * W.rhoE_ampa));
            IampaInp = single(full( ...
                W.S_Ilgn * lgnIEvents(:, randomBins(2,:)) + ...
                W.S_amb * ambIEvents(:, randomBins(4,:)) + ...
                W.S_IL6 * l6IEvents(:, randomBins(6,:)) * W.rhoI_ampa));
            EnmdaInp = single(full(W.S_EL6 * l6EEvents(:, randomBins(5,:)) * W.rhoE_nmda));
            InmdaInp = single(full(W.S_IL6 * l6IEvents(:, randomBins(6,:)) * W.rhoI_nmda));
        end

        [RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R, ...
            GE_ampa_D,GE_nmda_D,GE_gaba_D, ...
            RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R, ...
            GI_ampa_D,GI_nmda_D,GI_gaba_D, ...
            EEDlyRcd,IEDlyRcd,EIDlyRcd] = ...
            V1NetworkUpdate_Ver3_Drive_PoissonRead_EEnIEnEIDelay( ...
            RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R, ...
            GE_ampa_D,GE_nmda_D,GE_gaba_D, ...
            RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R, ...
            GI_ampa_D,GI_nmda_D,GI_gaba_D, ...
            EEDlyRcd,IEDlyRcd,EIDlyRcd, C_EE,C_EI,C_IE,C_II, ...
            W.S_EE,W.S_EI,W.S_IE,W.S_II, ...
            W.tau_ampa_R,W.tau_ampa_D,W.tau_nmda_R,W.tau_nmda_D, ...
            W.tau_gaba_R,W.tau_gaba_D,W.tau_ref,dt,W.p_EEFail, ...
            W.gL_E,W.Ve,W.rhoE_ampa,W.rhoE_nmda, ...
            W.gL_I,W.Vi,W.rhoI_ampa,W.rhoI_nmda, ...
            EampaInp(:,frameIndex),IampaInp(:,frameIndex), ...
            EnmdaInp(:,frameIndex),InmdaInp(:,frameIndex));

        globalTimeMs = ((section - 1) * stepsPerSection + timeN) * dt;
        if globalTimeMs >= analysisStartMs
            idxE = find(SpE); idxI = find(SpI);
            spikeCountE(idxE) = spikeCountE(idxE) + 1;
            spikeCountI(idxI) = spikeCountI(idxI) + 1;
        end
        if sum(isnan(VE)) > 0.8 * N_E
            error('Paper3 SNN became nonfinite at %.3f ms.', globalTimeMs);
        end
    end
    fprintf('angle %.2f replicate %d section %d/%d: %.1f s\n', ...
        angleDeg, replicate, section, TPar, toc(timerSection));
end

FrE = spikeCountE * (1000 / analysisMs);
FrI = spikeCountI * (1000 / analysisMs);
[FrSMap,NnSPixel] = NeuVec2Pixel(FrE(~EcplxInd), W.NnS, NPixX*N_HC, NPixY*N_HC);
[FrCMap,NnCPixel] = NeuVec2Pixel(FrE(EcplxInd), W.NnC, NPixX*N_HC, NPixY*N_HC);
[FrIMap,NnIPixel] = NeuVec2Pixel(FrI, W.NnI, NPixX*N_HC, NPixY*N_HC);

L4EE = C_EE * FrE; L4EI = C_EI * FrI;
L4IE = C_IE * FrE; L4II = C_II * FrI;
[L4SE,~] = NeuVec2Pixel(L4EE(~EcplxInd), W.NnS, NPixX*N_HC, NPixY*N_HC);
[L4CE,~] = NeuVec2Pixel(L4EE(EcplxInd), W.NnC, NPixX*N_HC, NPixY*N_HC);
[L4IEMap,~] = NeuVec2Pixel(L4IE, W.NnI, NPixX*N_HC, NPixY*N_HC);
[L4SI,~] = NeuVec2Pixel(L4EI(~EcplxInd), W.NnS, NPixX*N_HC, NPixY*N_HC);
[L4CI,~] = NeuVec2Pixel(L4EI(EcplxInd), W.NnC, NPixX*N_HC, NPixY*N_HC);
[L4IIMap,~] = NeuVec2Pixel(L4II, W.NnI, NPixX*N_HC, NPixY*N_HC);

pathwayScale = canonicalPathwayScales(W);
L4EByPopulation = [L4SE(:)/pathwayScale.L4SEp, ...
    L4CE(:)/pathwayScale.L4CEp, L4IEMap(:)/pathwayScale.L4IEp];
L4IByPopulation = [L4SI(:)/pathwayScale.L4SIp, ...
    L4CI(:)/pathwayScale.L4CIp, L4IIMap(:)/pathwayScale.L4IIp];

[L6SMap,~] = NeuVec2Pixel(rS_L6 * (1e3/W.NS_L6), W.NnS, NPixX*N_HC, NPixY*N_HC);
[L6CMap,~] = NeuVec2Pixel(rC_L6 * (1e3/W.NC_L6), W.NnC, NPixX*N_HC, NPixY*N_HC);
[L6IMap,~] = NeuVec2Pixel(rI_L6 * (1e3/W.NI_L6), W.NnI, NPixX*N_HC, NPixY*N_HC);
L6ByPopulation = [L6SMap(:), L6CMap(:), L6IMap(:)] / 3;

LGNWeightsS = pixelCategoryWeights(W.OD_S, NnSPixel, 4);
LGNWeightsC = pixelCategoryWeights(W.OD_C, NnCPixel, 4);
LGNWeightsI = pixelCategoryWeights(W.OD_I, NnIPixel, 4);

NWSmlt = struct();
NWSmlt.FS = FrSMap(:); NWSmlt.FC = FrCMap(:); NWSmlt.FI = FrIMap(:);
NWSmlt.L4EByPopulation = L4EByPopulation;
NWSmlt.L4IByPopulation = L4IByPopulation;
NWSmlt.L6ByPopulation = L6ByPopulation;
NWSmlt.LGNWeightsS = LGNWeightsS;
NWSmlt.LGNWeightsC = LGNWeightsC;
NWSmlt.LGNWeightsI = LGNWeightsI;

metadata = struct('angleDeg',angleDeg,'replicate',replicate,'seed',seed, ...
    'durationMs',durationMs,'analysisMs',analysisMs,'dtMs',dt, ...
    'runtimeSeconds',toc(timerAll),'networkE',N_E,'networkI',N_I, ...
    'connectionFile',connectionFile,'workspaceFile',workspaceFile, ...
    'model','Paper3_TestL6ShapeDriveHPC_Final', ...
    'L6Curve','L6CurveFinalize','L6OrdFHz',L6OrdFHz, ...
    'L6smear',L6smear,'L6trunc',L6trunc,'pathwayScale',pathwayScale);
outputFile = fullfile(outputDir, sprintf( ...
    'Paper3NWSimulationPix_%.1fdeg_rep%02d_h96_validation.mat', angleDeg, replicate));
save(outputFile, 'NWSmlt', 'metadata', '-v7.3');
fprintf('Saved %s\n', outputFile);
end

function scales = canonicalPathwayScales(W)
L4SE = W.C_SS_Pixel_Us*W.FrSPixVec + W.C_SC_Pixel_Us*W.FrCPixVec;
L4CE = W.C_CS_Pixel_Us*W.FrSPixVec + W.C_CC_Pixel_Us*W.FrCPixVec;
L4IE = W.C_IS_Pixel_Us*W.FrSPixVec + W.C_IC_Pixel_Us*W.FrCPixVec;
L4SI = W.C_SI_Pixel_Us*W.FrIPixVec;
L4CI = W.C_CI_Pixel_Us*W.FrIPixVec;
L4II = W.C_II_Pixel_Us*W.FrIPixVec;
L4Eall = L4SE + L4CE + L4IE;
L4Iall = L4SI + L4CI + L4II;
scales = struct('L4SEp',mean(L4SE./L4Eall),'L4CEp',mean(L4CE./L4Eall), ...
    'L4IEp',mean(L4IE./L4Eall),'L4SIp',mean(L4SI./L4Iall), ...
    'L4CIp',mean(L4CI./L4Iall),'L4IIp',mean(L4II./L4Iall));
end

function weights = pixelCategoryWeights(neuronCategory, pixelMap, categoryCount)
pixelCount = max(pixelMap.Vec);
weights = zeros(pixelCount, categoryCount);
for pixel = 1:pixelCount
    members = (pixelMap.Vec == pixel);
    for category = 1:categoryCount
        weights(pixel,category) = mean(neuronCategory(members) == category);
    end
end
assert(max(abs(sum(weights,2)-1)) < 1e-12, 'Pixel category weights do not sum to one.');
end

function value = envNumber(name, defaultValue)
textValue = getenv(name);
if isempty(textValue)
    value = defaultValue;
else
    value = str2double(textValue);
    if ~isfinite(value); error('Invalid numeric environment variable %s.', name); end
end
end
