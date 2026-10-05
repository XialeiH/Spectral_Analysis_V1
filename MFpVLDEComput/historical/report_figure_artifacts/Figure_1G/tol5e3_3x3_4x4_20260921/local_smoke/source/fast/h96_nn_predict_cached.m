function y = h96_nn_predict_cached(Xn, model)
% Evaluate one h96 population MLP from already-normalized features.
h = tanh(Xn * model.W1.' + model.b1);
h = tanh(h * model.W2.' + model.b2);
h = tanh(h * model.W3.' + model.b3);
z = h * model.W4.' + model.b4;
y = log1p(exp(-abs(z))) + max(z, 0);
end
