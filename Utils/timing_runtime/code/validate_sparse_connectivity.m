function validate_sparse_connectivity(runRoot)
dataFile = fullfile(runRoot,'repo','Data','Paper2_NetworkTuning','Fig1V4', ...
    'AllMFPixPara_Paper2TuneFig1V4D2.mat');
D = load(dataFile,'C_SS_Pixel_Us','N_HC');
legacyDir = fullfile(runRoot,'legacy_utils');
sparseDir = fullfile(runRoot,'large_utils');
addpath(legacyDir,'-begin');
legacy = AveSpatKer_Rec(D.C_SS_Pixel_Us,D.N_HC,4,4,10,10,1,0);
clear AveSpatKer_Rec
rmpath(legacyDir);
addpath(sparseDir,'-begin');
optimized = AveSpatKer_Rec(D.C_SS_Pixel_Us,D.N_HC,4,4,10,10,1,0);
difference = full(legacy)-full(optimized);
fprintf('CONNECTIVITY_MAX_ABS=%.17g\n',max(abs(difference),[],'all'));
fprintf('CONNECTIVITY_FRO_REL=%.17g\n',norm(difference,'fro')/norm(legacy,'fro'));
assert(max(abs(difference),[],'all')<1e-12,'Sparse connectivity mismatch.');
end
