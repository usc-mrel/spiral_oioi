function outputParams = final_wf_decomp_no_filter(algoParams,Nx,Ny,Nz,filtsize,datafull,BWfull,fm)

Ims_final = zeros(Nx,Ny,Nz,algoParams.M);
Ims_temp  = zeros(Nx,Ny,algoParams.M * algoParams.N,Nz);
for slicenum = 1:Nz
    if Nz == 1
        Ims = datafull;
    else
        Ims = datafull(:,:,slicenum,:);
    end
    for nn = 1:algoParams.N
        Ims(:,:,nn) = Ims(:,:,nn).* exp(-1i*2*pi*algoParams.te(nn).* fm(:,:,slicenum));
    end
    
    Ims_all(:,:,:,slicenum) = Ims;
    
    [rr,cc]=find(BWfull(:,:,slicenum));
    for k = 1:length(rr)
        S_hat = zeros(algoParams.N*2,1);
        S_hat(1:algoParams.N) = real(Ims(rr(k),cc(k),:));
        S_hat((algoParams.N+1):(2*algoParams.N)) = imag(Ims(rr(k),cc(k),:));
        p_hat = algoParams.Ainv * S_hat;
        for j = 1:algoParams.M
            Ims_final(rr(k),cc(k),slicenum,j) = p_hat((j*2-1)) + 1i*p_hat((j*2));
        end
        
        %% for error
        
        s_hat = algoParams.A * p_hat;
        Ims_temp(rr(k),cc(k),:,slicenum) = s_hat;
         
    end
end

outputParams.water = Ims_final(:,:,:,1);
outputParams.fat = Ims_final(:,:,:,2);
outputParams.fieldmap = fm;

outputParams.error = Ims_temp(:, :, 1:algoParams.N,:) + 1i * Ims_temp(:, :, algoParams.N + 1 : end, :);
outputParams.error = outputParams.error - Ims_all;

% outputParams.species(1).amps = Ims_final(:,:,:,1);
% outputParams.species(2).amps = Ims_final(:,:,:,2);