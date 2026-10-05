function eta = GA_nonnormality_index(A)
% eta = ||A'A - AA'||_F / ||A||_F^2 from proposal section 4.5.
denom = norm(A, 'fro')^2;
if denom == 0
    eta = 0;
    return
end
eta = norm(A' * A - A * A', 'fro') / denom;
end
