function x = GA_pack_state(state)
% Pack an LDE state struct into [S; C; I].
if isstruct(state)
    x = [state.S(:); state.C(:); state.I(:)];
else
    x = state(:);
end
end
