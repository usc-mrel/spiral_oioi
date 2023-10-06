function imDataParams = get_resolution_levels(imDataParams, l)

imDataParams.images_levels{1} = imDataParams.images;
if l == 1 
    return
end
for i = 2:l
    im_temp = imDataParams.images_levels{i-1};
    [sx, sy, sz, nc, ne] = size(im_temp);
    if mod(sx, 2) ~= 0
        im_temp(end+1, :, :, :, :) = 0;
        sx = sx + 1;
    end
    if mod(sy, 2) ~= 0
        im_temp(:, end+1, :, :, :) = 0;
        sy = sy + 1;
    end
    imDataParams.images_levels{i} = reshape(im_temp, [2, sx/2, 2, sy/2, sz, nc, ne]);
    imDataParams.images_levels{i} = sum(imDataParams.images_levels{i}, 1);
    imDataParams.images_levels{i} = sum(imDataParams.images_levels{i}, 3);
    imDataParams.images_levels{i} = reshape(imDataParams.images_levels{i}, [sx/2, sy/2, sz, nc, ne]);
end
