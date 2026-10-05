function replot_with_fixed_point_row(dataFile,setupFile,outputPdf,enrichedDataFile)
% Append an exact fixed-point snapshot and replot the E-only perturbation page.

arguments
    dataFile (1,:) char
    setupFile (1,:) char
    outputPdf (1,:) char
    enrichedDataFile (1,:) char
end

loadedResult = load(dataFile,'output');
output = loadedResult.output;
loadedSetup = load(setupFile,'setup');
setup = loadedSetup.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);

jFixed = l6ns_state_jacobian(setup,'L6',fixed,output.L6Weight);
[returnMode,~] = eigs(jFixed,1,'largestreal');
returnMode = local_orient_vector(returnMode,n,context.CWeight);
[singularOutput,singularValue,singularInput] = svds(jFixed,1,'largest');
singularInput = local_orient_vector(singularInput,n,context.CWeight);
if real(singularOutput'*singularInput)<0
    singularOutput = -singularOutput;
end

for index = 1:numel(output.AllSnapshots)
    output.AllSnapshots(index).PhaseName = '';
end
snapshot = output.AllSnapshots(end);
snapshot.TimeMs = output.Times(end);
snapshot.State = fixed;
snapshot.HCnorm = 0;
snapshot.DeltaNorm = 0;
snapshot.SingularAlignment = 0;
snapshot.ReturnAlignment = 1;
snapshot.JacobianMaxReal = real(output.FixedPointLeadingEigenvalue);
snapshot.Eigenvalues = output.FixedPointEigenvalues;
snapshot.ReturnDirectionMaps = local_signed_maps(returnMode,mapSize,context.CWeight);
snapshot.SingularInputMaps = local_signed_maps(singularInput,mapSize,context.CWeight);
snapshot.SingularOutputMaps = local_signed_maps(singularOutput,mapSize,context.CWeight);
snapshot.FiringMaps = local_state_maps(fixed,mapSize,context.CWeight);
snapshot.TopSingularValue = singularValue;
snapshot.PhaseName = 'Fixed point';
output.AllSnapshots(end+1) = snapshot;

save(enrichedDataFile,'output','-v7.3');
replot_perturbation_e_only(enrichedDataFile,outputPdf);
end

function maps = local_state_maps(state,mapSize,wC)
n = prod(mapSize);
s = reshape(state(1:n),mapSize);
c = reshape(state(n+(1:n)),mapSize);
i = reshape(state(2*n+(1:n)),mapSize);
e = (1-wC)*s+wC*c;
maps = {s;c;i;e};
end

function maps = local_signed_maps(vector,mapSize,wC)
maps = local_state_maps(real(vector),mapSize,wC);
for index = 1:numel(maps)
    maps{index} = maps{index}/max(norm(maps{index}(:),2),eps);
end
end

function vector = local_orient_vector(vector,n,wC)
e = (1-wC)*real(vector(1:n))+wC*real(vector(n+(1:n)));
[~,index] = max(abs(e));
if e(index)<0
    vector = -vector;
end
end
