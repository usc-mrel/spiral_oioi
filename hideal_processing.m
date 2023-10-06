ccc

format compact;
warning ('off','all');

%% ADD PATH
addpath ./mfile/hu/
addpath ./mfile/hu/'Support Functions'/
% BASEPATH = '/server/sdata/ytian/water_fat/fwtoolbox_v1_code/hu/';
% tmp = BASEPATH; addpath(tmp);
% tmp = fullfile(BASEPATH,'Support Functions'); addpath(tmp);
% tmp = fullfile(BASEPATH,'Data'); addpath(tmp);
% tmp = fullfile(BASEPATH,'..', '..', 'fwtoolbox_v1_data', 'USC'); addpath(tmp);
% tmp = fullfile(tmp,'SingleChannel 3 Echo');addpath(tmp)

%%
all_mat = dir('./recon_data/*out_in_out_in*.mat');
nfile = length(all_mat);

for ifile = 1:nfile

load(fullfile(all_mat(ifile).folder, all_mat(ifile).name), 'im_king', 'para')
tic
TE = para.kspace_info.user_TE / 1000;
EchoSpacing = para.kspace_info.user_readoutTime / 2 / 1000;

im_king = im_king(:, :, 11:end-10 , :);
imDataParams.images = mean(im_king, 3);

imDataParams.images = permute(imDataParams.images, [1, 2, 3, 5, 4]);
imDataParams.TE = [TE, TE + EchoSpacing, TE + EchoSpacing, TE + EchoSpacing * 2];
imDataParams.FieldStrength = 0.55;
imDataParams.PrecessionIsClockwise = 0;

clear image*

%% time direction
if imDataParams.PrecessionIsClockwise <= 0
    imDataParams.images = conj(imDataParams.images);
end

%% parameters 
algoParams.M = 2;% Number of species; 2 = fat/water

% water
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;

% fat
algoParams.species(2).name = 'fat';
fat_ppm = -[3.8 3.4 2.6 1.94 0.39 -0.6]; % negative fat frequency --> 2*pi*fieldmap*te effect
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.species(2).frequency = fat_ppm .* 4258 * imDataParams.FieldStrength * 1e-2;
algoParams.N = length(imDataParams.TE);
algoParams.te = imDataParams.TE;
[algoParams.A, algoParams.C, algoParams.D] = solveA(algoParams);
algoParams.Ainv = (transpose(algoParams.A)*algoParams.A) \ transpose(algoParams.A);

algoParams.maxiter      = 300;
algoParams.fm_epsilon   = 1;

%% resolution levels
l = 5;
imDataParams = get_resolution_levels(imDataParams, l);

BW = sos(imDataParams.images_levels{l}) > 0;
[sx, sy, sz, nc, ne] = size(imDataParams.images_levels{l});

%% initial guess of field map
fm_hold = run_ideal(algoParams, sx, sy, sz, BW, squeeze(imDataParams.images_levels{l}));

%% regional grow
imDataParams.fm_levels{l} = regional_grow_circular(squeeze(imDataParams.images_levels{l}), algoParams, fm_hold);

%% iterate layers
for i = l-1:-1:1
    sx = size(imDataParams.images_levels{i}, 1);
    sy = size(imDataParams.images_levels{i}, 2);
    fm_hold = imresize(imDataParams.fm_levels{i+1}, [sx, sy]);
    BW = sos(imDataParams.images_levels{i}) > 0;
    fm_hold = run_ideal(algoParams, sx, sy, 1, BW, squeeze(imDataParams.images_levels{i}), fm_hold);
    fm_hold = medfilt2(fm_hold, [max(round(sx * 0.02), 3), max(round(sy * 0.02), 3)]);
    w = get_b0_filter(sx, sy);
    fm_hold = real(ifft2(fftshift2(fftshift2(fft2(fm_hold)) .* w)));
    imDataParams.fm_levels{i} = fm_hold;
    imageMRI(fm_hold);
end
outputParams_hideal_averge = final_wf_decomp_no_filter(algoParams, sx, sy, sz, 1, squeeze(imDataParams.images_levels{1}), BW, fm_hold);
f = imageMRI([outputParams_hideal_averge.water, outputParams_hideal_averge.fat]);
brighten(0.2)
% exportgraphics(f, ['./figures/', all_mat(ifile).name(1:end-4), '_king_hideal.png'], 'Resolution', 300)
cc

%% apply on all frames
nframe = size(im_king, 3);
BW = sos(im_king) > max(vec(sos(im_king))) * 0.02;
fm_all = run_ideal(algoParams, sx, sy, nframe, BW, conj(im_king), repmat(fm_hold, [1, 1, nframe]));

fm_smooth = medfilt3(fm_all, [7,7,7]);
fm_smooth = real(ifft2(fftshift2(fftshift2(fft2(fm_smooth)) .* w)));
outputParams_hideal = final_wf_decomp_no_filter(algoParams, sx, sy, nframe, 1, conj(im_king), BW, fm_smooth);

outputParams_hideal.fm_estimate     = single(fm_all);
outputParams_hideal.water           = single(outputParams_hideal.water);
outputParams_hideal.fat             = single(outputParams_hideal.fat);
outputParams_hideal.fieldmap        = single(outputParams_hideal.fieldmap);
outputParams_hideal.error           = single(outputParams_hideal.error);
toc
if ~isfolder('./hideal/')
    mkdir hideal
end
save(fullfile('./hideal/', all_mat(ifile).name), 'outputParams_hideal', 'outputParams_hideal_averge', '-v7.3')
end



