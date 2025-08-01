function fm_est = regional_grow_circular(im, algoParams, fm_hold)
sx = size(im, 1);
sy = size(im, 2);
r1 = 1;
r2 = 2;

start_pt = round([sx/2, sy/2]);

mask = false(sx, sy);
mask(start_pt(1), start_pt(2)) = true;
mask = round(bwdist(mask));
fm_est = zeros(sx, sy);
mask_ = mask<=r1;
fm_est(mask_) = fm_hold(mask_);

mask_run = false(sx, sy);
mask_run(mask_) = true;

for ii = r1+1:max(mask(:))
    mask_ = mask == ii;
    [xx, yy] = find(mask_);
    n = length(xx);
    fm_temp = zeros(sx, sy);
    for jj = 1:n
        %         if mask_im(xx(jj), yy(jj))
        %             if xx(jj) == 27  && yy(jj) == 23
        %                 keyboard
        %             end
        pixel = im(xx(jj), yy(jj), :);
        
        for kk = -123:123
            outputParams = final_wf_decomp_no_filter(algoParams, 1, 1, 1, 1, pixel, 1, kk);
            error(kk + 124) = sos(outputParams.error);
        end
        local_mins = local_max(-error) - 124;
        local_mins(local_mins == 123) = [];
        local_mins(local_mins == -123) = [];
        if isempty(local_mins)
            local_mins = [-123, 123];
        end
        
        fm_candidate = [local_mins - 247, local_mins, local_mins + 247];
        
        x_0 = max(1,  xx(jj) - r2);
        x_1 = min(sx, xx(jj) + r2);
        y_0 = max(1,  yy(jj) - r2);
        y_1 = min(sy, yy(jj) + r2);
        
        box = sos(im(x_0:x_1, y_0:y_1, :)) .* mask_run(x_0:x_1, y_0:y_1);
        box = box ./ sum(sum(box));
        
        fm_ave = sum(sum(box .* fm_est(x_0:x_1, y_0:y_1)));
        
        [~, idx] = min(abs(fm_ave - fm_candidate));
        
        fm_est(xx(jj), yy(jj)) = fm_candidate(idx);
        mask_run(xx(jj), yy(jj)) = true;
        %         end
    end
    %     for jj = 1:n
    %         if mask_im(xx(jj), yy(jj))
    %             mask_run(xx(jj), yy(jj)) = true;
    %         end
    %     end
    %     fm_est = fm_est + fm_temp;
    
end
% imageMRI(fm_est);
% drawnow
