function state = GA_unpack_state(x, template)
% Unpack [S; C; I] using an existing LDE state as shape template.
n = numel(template.S);
state = template;
state.S = reshape(x(1:n), size(template.S));
state.C = reshape(x(n+(1:n)), size(template.C));
state.I = reshape(x(2*n+(1:n)), size(template.I));
end
