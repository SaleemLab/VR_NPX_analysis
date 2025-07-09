%%%%%%% Main Time frequency analysis (Average activities + extract -200-0 and 0-200ms phase coupling and average SO/Spindle/ripples power)

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))



%% Visualising TF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(analysis_folder,'ripples_TF_stats_POST.mat'),'ripples_TF_stats')
load(fullfile(analysis_folder,'UP_TF_stats_POST.mat'),'UP_TF_stats')
load(fullfile(analysis_folder,'DOWN_TF_stats_POST.mat'),'DOWN_TF_stats')
load(fullfile(analysis_folder,'spindles_TF_stats_POST.mat'),'spindles_TF_stats')

timevec = ripples_TF_stats.timebin{1};
freqs  = ripples_TF_stats.freq{1};

fig = figure('Name', ['TF amp ripples all POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,ripples_TF_stats.mean_V1_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp Ripples all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp Ripples all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_HPC_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)






fig = figure('Name', ['TF amp ripples DOWN POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,ripples_TF_stats.mean_V1_amp_DOWN{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp Ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp Ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_HPC_amp_DOWN{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)





fig = figure('Name', ['TF amp ripples UP POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,ripples_TF_stats.mean_V1_amp_UP{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp Ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp Ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_HPC_amp_UP{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp Ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-6 6])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)






fig = figure('Name', ['TF amp UP POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,UP_TF_stats.mean_V1_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,UP_TF_stats.mean_HPC_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)





fig = figure('Name', ['TF amp UP ripples POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,UP_TF_stats.mean_V1_amp_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,UP_TF_stats.mean_HPC_amp_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp UP ripple')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)





fig = figure('Name', ['TF amp DOWN POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,DOWN_TF_stats.mean_V1_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,DOWN_TF_stats.mean_HPC_amp_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)







fig = figure('Name', ['TF amp DOWN ripples POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(4,DOWN_TF_stats.mean_V1_amp_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));
[tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 amp DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 amp DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,DOWN_TF_stats.mean_HPC_amp_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi HPC amp DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
xlim([-1.5 1.5])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)

save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','TF analysis'),[],'ContentType','vector')







fig = figure('Name', ['TF PLV DOWN all POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,DOWN_TF_stats.mean_PLV_V1_all{:});
ipsi_map = tf_map;
tf_map = cat(3,DOWN_TF_stats.mean_PLV_HPC_all{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,DOWN_TF_stats.mean_PLV_V1_HPC_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV DOWN all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV DOWN all')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)









fig = figure('Name', ['TF PLV DOWN ripples POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,DOWN_TF_stats.mean_PLV_V1_ripples{:});
ipsi_map = tf_map;
tf_map = cat(3,DOWN_TF_stats.mean_PLV_HPC_ripples{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,DOWN_TF_stats.mean_PLV_V1_HPC_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV DOWN ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)







fig = figure('Name', ['TF PLV UP all POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,UP_TF_stats.mean_PLV_V1_all{:});
ipsi_map = tf_map;
tf_map = cat(3,UP_TF_stats.mean_PLV_HPC_all{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,UP_TF_stats.mean_PLV_V1_HPC_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)






fig = figure('Name', ['TF PLV UP ripples POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,UP_TF_stats.mean_PLV_V1_ripples{:});
ipsi_map = tf_map;
tf_map = cat(3,UP_TF_stats.mean_PLV_HPC_ripples{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,UP_TF_stats.mean_PLV_V1_HPC_ripples{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV UP ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)




fig = figure('Name', ['TF PLV ripples all POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,ripples_TF_stats.mean_PLV_V1_all{:});
ipsi_map = tf_map;
tf_map = cat(3,ripples_TF_stats.mean_PLV_HPC_all{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_PLV_V1_HPC_all{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV ripples')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)











fig = figure('Name', ['TF PLV ripples UP POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,ripples_TF_stats.mean_PLV_V1_UP{:});
ipsi_map = tf_map;
tf_map = cat(3,ripples_TF_stats.mean_PLV_HPC_UP{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_PLV_V1_HPC_UP{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV ripples UP')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)







fig = figure('Name', ['TF PLV ripples UP POST'], ...
    'Position', [100 100 1200 500]);
tf_map = cat(3,ripples_TF_stats.mean_PLV_V1_DOWN{:});
ipsi_map = tf_map;
tf_map = cat(3,ripples_TF_stats.mean_PLV_HPC_DOWN{:});
contra_map =  tf_map;
% [tf_mean, sig_mask] = permutation_TF_test(ipsi_map, 1000);

subplot(2,2,1)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('V1 PLV ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
% clim([-7 3])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)
%             sig_mask = curr_TF_stats.(['mean_PLV_sig_mask_V1', suffix]){nsession};
% hold on
% contour(timevec, log2(freqs), sig_mask, [1 1], 'r', 'LineWidth', 1.2)

subplot(2,2,2)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('HPC PLV ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
% clim([-3.14 3.14])
ylim([min(log2(freqs)) max(log2(freqs))])
clim([0 0.5])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



tf_map = cat(4,ripples_TF_stats.mean_PLV_V1_HPC_DOWN{:});
ipsi_map = squeeze(tf_map(1,:,:,:));
contra_map =  squeeze(tf_map(2,:,:,:));

subplot(2,2,3)
tf_mean = mean(ipsi_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('ipsi V1 HPC PLV ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)


subplot(2,2,4)
tf_mean = mean(contra_map,3,'omitnan');
%             contourf(timevec, log2(freqs), tf_mean, 40, 'linecolor', 'none')
imagesc(timevec, log2(freqs), tf_mean)
axis xy; title('contra V1 HPC PLV ripples DOWN')
xlabel('Time (s)'), ylabel('log_2(Freq Hz)')
set(gca, 'YTick', log2([1 2 4 8 16 32 64 128 256]), ...
    'YTickLabel', {'1','2','4','8','16','32','64','128','256'})
colorbar
xline(0,'r--')
clim([0 0.5])
ylim([min(log2(freqs)) max(log2(freqs))])
yline(log2([9 17])); yline(log2([1 4])); yline(log2([125 300]))
set(gca, 'TickDir','out', 'Box','off', 'FontSize',12)



save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','TF analysis'),[],'ContentType','vector')
