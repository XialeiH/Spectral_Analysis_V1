function d2 = GA_local_discriminability(deltaCoord, metric)
% Local squared distinguishability ds^2 = deltaCoord' * g * deltaCoord.
deltaCoord = deltaCoord(:);
d2 = deltaCoord' * metric * deltaCoord;
end
