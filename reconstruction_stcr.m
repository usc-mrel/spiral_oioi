ccc
addpath ./mfile/

all_mat = dir('./raw/*.mat');
nfile = length(all_mat);

%% set settings
para.setting.ifplot     = 0;                    % plot convergence during reconstruction 
para.setting.ifGPU      = 1;                    % set to 1 when you want to use GPU

%% set recon parameters
para.weight_tTV         = 0.01;                 % temporal TV regularizaiton parameter (normalized by F^T d)
para.weight_sTV         = 0.001;                % spatial TV regularizaiton parameter (normalized by F^T d)

para.Recon.epsilon      = eps('single');        % small vale to avoid singularity in TV constraint
para.Recon.step_size    = 2;                    % initial step size
para.Recon.noi          = 10;                   % number of NLsCG iterations
para.Recon.type         = '2D Spiral server';   % 2D spiral
para.Recon.break        = 1;                    % stop iteration if creteria met. Otherwise will run to noi

para.Recon.FOV_recon            = [560, 560];   % mm
para.Recon.temporal_resolution  = 45;           % ms

%% recon
for i = 1:nfile
    para.raw_file  = fullfile(all_mat(i).folder, all_mat(i).name);
    para.file_name = all_mat(i).name;
    
    recon_spiral_oioi(para)
end