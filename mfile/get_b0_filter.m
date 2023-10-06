function weights = get_b0_filter(Nx, Ny)

w_fB0 = 64;  % 32 in BART?
h_fB0 = 8;   % Sobolev index for B0 field inhomogeneity
weights = zeros(Nx, Ny, 'single');
for idx2 = 1:Nx
    for idx1 = 1:Ny
        %------------------------------------------------------------------
        % Calculate the k-space weight for B0 field inhomogeneity
        %------------------------------------------------------------------
        kx = (-floor(Nx/2) + idx1 - 1) / Nx;
        ky = (-floor(Ny/2) + idx2 - 1) / Ny;
        weights(idx1,idx2) = 1 / (1 + w_fB0 * (kx^2 + ky^2))^h_fB0;
    end
end