clear;

file = "F:\TEAMI\Colum_K3_hBN\pos_refine\MLs_p8_g64_nv128_betaO0.1_betaP0.1_mm_Ns17_dz10_reg0.8_fixedPos_1\Niter200.mat";
load(file);
% load("D:\0_Codes\ePIE-Code_Miao\for_Colum_Haozhi\annealing_pos_correction_jianhua\offset_probe_true.mat")

% pos_t = squeeze(cat(3, offset_probe_b_true, offset_probe_a_true));
pos0 = outputs.probe_positions_0;
pos = outputs.probe_positions;
% cent_t = mean(pos_t, 1);
cent0 = mean(pos0, 1);
cent = mean(pos, 1);
fprintf('center of pos0: %.3e, %.3e\ncenter of pos: %.3e, %.3e\n', cent0, cent)
% pos_t = pos_t + cent0 - cent_t;
% pos = pos + cent0 - cent;

f = figure(1);
f.Position = [500 500 800 800];

% subplot(131)
% scatter(pos_t(:,1), pos_t(:,2), 40, 'black'); hold on;
s = scatter(pos(:,1), pos(:,2), 10, 'black', "filled");hold on;
% s.AlphaData = 0.5*ones(size(pos, 1), 1);
% s.MarkerFaceAlpha = 'flat';
scatter(pos0(:,1), pos0(:,2), 5, [0.9290 0.6940 0.1250], 'filled'); hold on;

% quiver(pos0(:,1), pos0(:,2), pos(:,1)-pos0(:,1), pos(:,2)-pos0(:,2), 'black')
legend('pos after opt.','input pos')
axis square

% dev = pos - pos0;
% subplot(132)
% scatter(1:size(dev, 1), dev(:,1), 'filled');
% % ylim([-0.6, 0.6])
% title('deviation in x');
% subplot(133)
% scatter(1:size(dev, 1), dev(:,2), 'filled');
% % ylim([-0.6, 0.6])
% title('deviation in y');