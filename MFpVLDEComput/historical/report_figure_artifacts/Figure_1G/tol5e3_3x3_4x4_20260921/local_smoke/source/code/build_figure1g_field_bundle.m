function build_figure1g_field_bundle(runRoot, fieldHC)
% Precompute all field-size-dependent connectivity outside timed solves.
fieldHC = str2double(string(fieldHC));
repoRoot = fullfile(runRoot, 'repo');
addpath(fullfile(repoRoot, 'Utils'), '-begin');
addpath(fullfile(runRoot, 'large_utils'), '-begin');
addpath(fullfile(runRoot, 'fast'), '-begin');

dataFile = fullfile(repoRoot, 'Data', 'Paper2_NetworkTuning', 'Fig1V4', ...
    'AllMFPixPara_Paper2TuneFig1V4D2.mat');
D = load(dataFile);

nPixX = 10;
nPixY = 10;
n = fieldHC^2 * nPixX * nPixY;
paraE = 1;
paraI = 1;
boundaryMode = 'periodic';

frS = D.FrSPixVec;
frC = D.FrCPixVec;
frI = D.FrIPixVec;
l4SE = D.C_SS_Pixel_Us * frS + D.C_SC_Pixel_Us * frC;
l4SI = D.C_SI_Pixel_Us * frI;
l4CE = D.C_CS_Pixel_Us * frS + D.C_CC_Pixel_Us * frC;
l4CI = D.C_CI_Pixel_Us * frI;
l4IE = D.C_IS_Pixel_Us * frS + D.C_IC_Pixel_Us * frC;
l4II = D.C_II_Pixel_Us * frI;
l4Eall = l4SE + l4CE + l4IE;
l4Iall = l4SI + l4CI + l4II;

ctx.L4SEp = mean(l4SE ./ l4Eall);
ctx.L4CEp = mean(l4CE ./ l4Eall);
ctx.L4IEp = mean(l4IE ./ l4Eall);
ctx.L4SIp = mean(l4SI ./ l4Iall);
ctx.L4CIp = mean(l4CI ./ l4Iall);
ctx.L4IIp = mean(l4II ./ l4Iall);

fprintf('Building %dx%d HC connectivity (%d pixels)\n', fieldHC, fieldHC, n);
ctx.C_SS_meanU = make_connection(D.C_SS_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_CS_meanU = make_connection(D.C_CS_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_IS_mean  = make_connection(D.C_IS_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_SC_meanU = make_connection(D.C_SC_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_CC_meanU = make_connection(D.C_CC_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_IC_mean  = make_connection(D.C_IC_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraE);
ctx.C_SI_mean  = make_connection(D.C_SI_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraI);
ctx.C_CI_mean  = make_connection(D.C_CI_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraI);
ctx.C_II_mean  = make_connection(D.C_II_Pixel_Us, D.N_HC, fieldHC, nPixX, nPixY, paraI);

cplxR = D.N_C / D.N_E;
ctx.C_SS_meanU = ctx.C_SS_meanU * ((D.N_SS - 6*(1-cplxR)) / D.N_SS);
ctx.C_SC_meanU = ctx.C_SC_meanU * ((D.N_SC - 4*cplxR) / D.N_SC);
ctx.C_CS_meanU = ctx.C_CS_meanU * ((D.N_CS - 6*(1-cplxR)) / D.N_CS);
ctx.C_CC_meanU = ctx.C_CC_meanU * ((D.N_CC - 4*cplxR) / D.N_CC);

lgnFiltered = SpatialGaussianFilt_my(D.OD_SMap, 3, D.n_S_HC, D.n_S_HC*0.2, 1, false);
lgnFiltered = [lgnFiltered, zeros(size(lgnFiltered,1),1)];
ctx.PixLGNCtgr = LGNIndSpat_Rec(lgnFiltered, 1:5, D.NnSPixel, D.N_HC, ...
    fieldHC, fieldHC, nPixX, nPixY, true);
ctx.PixLGNCtgr = Ocu_LGNL6symm(ctx.PixLGNCtgr, fieldHC, fieldHC, nPixX, nPixY);

ctx.N_HCOutY = fieldHC;
ctx.NPixX = nPixX;
ctx.NPixY = nPixY;
ctx.p = 0.33;
ctx.L6Kernel = make_l6_kernel();
ctx.L6parId = 5;
ctx.L6pars = cell(1,5);
ctx.L6pars{5} = {2.5,3,62.5,86,80.4,106,90,142.2,102,166.7,108,210,118, ...
    {[0 15 29],[0,34.0],[28 32 50],[69 85]}};
ctx.EKpUse = {1,0,[50;70;100],[50;60.5;69],'quadratic'};
ctx.IKpUse = struct('Thrsld1',54,'Thrsld2',63,'Highist',103,'Slope',0.941, ...
    'down1',0.3,'IntaH',4.5,'IntaL',-2.5,'IntbL',4,'IntbH',-0.3, ...
    'Mode','multisigmoid');

initialState = struct('S', 2.5*ones(n,1), 'C', 8*ones(n,1), 'I', 18*ones(n,1));
bundleFile = fullfile(runRoot, 'bundles', sprintf('field_%02dHC.mat', fieldHC));
save(bundleFile, 'ctx', 'initialState', 'fieldHC', 'boundaryMode', '-v7.3');
fprintf('Saved %s\n', bundleFile);
end

function A = make_connection(source, nHCIn, nHCOut, nPixX, nPixY, para)
A = sparse(AveSpatKer_Rec(source, nHCIn, nHCOut, nHCOut, ...
    nPixX, nPixY, para, 0));
end

function kernel = make_l6_kernel()
sigma = 0.75;
truncation = 1.5 / sigma;
kernelSize = max(2, floor(2*sigma*truncation));
[x,y] = meshgrid(1:kernelSize, 1:kernelSize);
center = kernelSize/2;
kernel = exp(((x-center).^2 + (y-center).^2) ./ (-2*sigma^2));
kernel(kernel < exp(-(truncation^2)/2)) = 0;
kernel = kernel + kernel(end:-1:1,:) + kernel(:,end:-1:1) + kernel(end:-1:1,end:-1:1);
kernel = kernel / sum(kernel,'all');
kernel = kernel / ((sum(kernel,'all')-kernel(2,2))/(1-0.3));
kernel(2,2) = 0.3;
end
