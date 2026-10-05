function spec = GA_h96_parameter_spec()
% Parameter directions for local h96 circuit-strength geometry.
% Names use target_from_source convention for the lateral matrices.
names = {'S_from_S','S_from_C','C_from_S','C_from_C', ...
    'I_from_S','I_from_C','S_from_I','C_from_I','I_from_I', ...
    'L6_feedback'};
groups = {'E_to_E','E_to_E','E_to_E','E_to_E', ...
    'E_to_I','E_to_I','I_to_E','I_to_E','I_to_I', ...
    'L6_to_L4'};

spec = struct();
spec.Names = names;
spec.Groups = groups;
spec.NumParameters = numel(names);
spec.Description = ['Local log-strength P directions. Lateral directions scale ', ...
    'one connectivity matrix; L6_feedback scales the L6 convolution input.'];
end
