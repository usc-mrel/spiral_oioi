ccc

all_mat = dir('./recon_data/*out_in_out_in*.mat');
nfile = length(all_mat);

for ifile = 1:nfile
   %% read data 
    file_name = fullfile(all_mat(ifile).folder, all_mat(ifile).name);
    variableInfo = who('-file', file_name);
    if ~ismember('im_king', variableInfo)
    load(file_name, 'image_out_1', 'image_in_1', 'image_out_2', 'image_in_2', 'para')
    t1 = tic;
    %% some informaiton
    kspace_info = para.kspace_info;
    res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
    
    [sx, sy, nframe] = size(image_out_1);
    matrix_size = [sx, sy];
    
    kx = kspace_info.kx_GIRF * matrix_size(1);
    ky = kspace_info.ky_GIRF * matrix_size(2);

    %% Calculate the time courses of the phase coefficients
    n_kappa = kspace_info.user_samples;
    
    ni = size(kx, 2);
    B0 = 0.55; % main field strength [T]
    
    % scale kx, ky to physical unit
    kxkymax = max(vec(sqrt(kx.^2 + ky.^2)));
    kxkymax_physical = 1 / res(1) * 1e3 * pi;
    
    kx_ = kx / kxkymax * kxkymax_physical; % [cycle/cm]
    ky_ = ky / kxkymax * kxkymax_physical; % [cycle/cm]
    
    gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]
    
    tdwell = 1 / kspace_info.user_samplingRate / 1000; % [sec]
    
    kx_ = kx_ / gamma / tdwell * 100;
    ky_ = ky_ / gamma / tdwell * 100;
    
    gz = zeros(n_kappa, ni, 'single'); % [G/cm]
    gx = [kx_(1, :); diff(kx_)];
    gy = [ky_(1, :); diff(ky_)];
    
    %% king's method
    opxres = matrix_size(1);
    opyres = matrix_size(2);
    
    dx = kspace_info.user_ResolutionX / 10; % [cm]
    dy = kspace_info.user_ResolutionY / 10; % [cm]
    
    rot = Q2Rot(kspace_info.user_QuaternionW, kspace_info.user_QuaternionX, kspace_info.user_QuaternionY, kspace_info.user_QuaternionZ);
    
    field_of_view_mm = matrix_size(1) * res(1);

    % this is the physical offset
    DCS_offset = [kspace_info.user_TranslationX; kspace_info.user_TranslationY; kspace_info.user_TranslationZ];
    % translate to logical offset
    DCS_offset = rot' * DCS_offset;
    
    % logical coordinate 
    x_range = (-floor(opxres/2):ceil(opxres/2)-1).' * dx; % [cm]
    y_range = (-floor(opyres/2):ceil(opyres/2)-1).' * dy;
    
    x = repmat(x_range', [matrix_size(2), 1]) / 100; % [m]
    y = repmat(y_range,  [1, matrix_size(1)]) / 100;
    z = zeros(matrix_size);

    x = x + DCS_offset(1) / 1000;
    y = y + DCS_offset(2) / 1000;
    z = z + DCS_offset(3) / 1000;
    
    image = cat(4, image_out_1, image_in_1, image_out_2, image_in_2);
    [im_king, Image_recon] = perform_deblurring_king_method_oioi(image, kx / matrix_size(1), ky / matrix_size(2), gx, gy, gz, x, y, z, field_of_view_mm, DCS_offset / 1e3, gamma, B0, tdwell, rot);
    time_cf = toc(t1);
    
    %% save
    save(fullfile(all_mat(ifile).folder, all_mat(ifile).name), 'im_king', 'time_cf', '-append')
    end
end
