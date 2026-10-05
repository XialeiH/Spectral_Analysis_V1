function K = build_conv_matrix_circular(kernel, Ny, Nx)
    Ny_in = double(Ny);
    Nx_in = double(Nx);
    Ny = round(Ny_in);
    Nx = round(Nx_in);
    tol = 1e-8;
    if ~isfinite(Ny_in) || ~isfinite(Nx_in) || Ny < 1 || Nx < 1
        error('build_conv_matrix_circular:InvalidSize', ...
            'Ny/Nx must be finite positive values. Got Ny=%.12g, Nx=%.12g.', Ny_in, Nx_in);
    end
    if abs(Ny_in - Ny) > tol || abs(Nx_in - Nx) > tol
        error('build_conv_matrix_circular:NonIntegerSize', ...
            ['Ny/Nx must be integer-valued (or within %.1e tolerance). ', ...
             'Got Ny=%.12g, Nx=%.12g.'], tol, Ny_in, Nx_in);
    end

    persistent cache;
    key = sprintf('k%dx%d_%dx%d', Ny, Nx, size(kernel,1), size(kernel,2));
    if ~isempty(cache) && isfield(cache, key)
        K = cache.(key);
        return;
    end

    N = Ny * Nx;
    idx_grid = reshape(1:N, Ny, Nx);

    [kh, kw] = size(kernel);
    ctr_r = ceil(kh/2);
    ctr_c = ceil(kw/2);

    rows = [];
    cols = [];
    vals = [];

    for r = 1:kh
        for c = 1:kw
            weight = kernel(r,c);
            if weight == 0, continue; end
            dr = r - ctr_r;
            dc = c - ctr_c;
            % conv2-style (flipped kernel) with circular padding
            src = circshift(idx_grid, [dr, dc]);
            rows = [rows; idx_grid(:)];
            cols = [cols; src(:)];
            vals = [vals; weight * ones(N,1)];
        end
    end

    K = sparse(rows, cols, vals, N, N);
    if isempty(cache)
        cache = struct();
    end
    cache.(key) = K;
end
