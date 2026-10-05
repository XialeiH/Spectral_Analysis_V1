function fast = h96_nn_set_condition_staticfirst(base,contrast,orientation)
% Precompute the static orientation/contrast contribution to layer 1.
fast = h96_nn_set_condition(base,contrast,orientation);
staticInputs = fast.staticCompactPages(:,1:2,:);
fast.staticLayer1 = pagemtimes(staticInputs,fast.W1T(1:2,:,:))+fast.b1;
fast.dynamicW1T = fast.W1T(3:5,:,:);
fast = rmfield(fast,{'staticCompactPages','W1T'});
end
