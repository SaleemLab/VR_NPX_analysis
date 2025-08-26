
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\fig2svg'))

addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


close all
hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\Ripple power predicted by spindle power.fig', 'visible');
figure(hFig)

% exportgraphics(gcf, 'Ripple power predicted by spindle power.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'none');

plot2svg('Ripple power predicted by spindle power rasterised.svg');

options.subplot_idx = [1 2];
smart_rasterise_scatter(hFig, 'Ripple power predicted by spindle power rasterised.pdf', struct('subplot_idx', [1 2]));'



hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\Ripple power predicted by spindle power.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\Ripple power predicted by spindle phase.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\Ripple power predicted by SO phase.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\Ripple power predicted by SO power.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (full windows)\Ripples predicted by spindles.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (full windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


% hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows)\UP-DOWN lag vs ripple lag.fig', 'visible');
% print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (full windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
% close all


% hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows)\UP DOWN lag predicted by ripples.fig', 'visible');
% print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (150ms windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
% close all


hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows)\Ripples predicted by DOWN UP synchrony.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (150ms windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows)\HPC excitation predicted by DOWN UP synchrony.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (150ms windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (full windows)\V1 depression predicted by ripples.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (full windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction'),[],'ContentType','vector')



hFig = openfig('P:\PhD\Dissertation\Figure materials\M24064 20241206 Checkerboard_sh2_half event (all filtered) probe 1 X shank 2 ADJUSTED', 'visible');

print(hFig, sprintf('P:/PhD/Dissertation/Figure materials/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows)\UP DOWN lag predicted by ripples.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (150ms windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (150ms windows UP lag control)\UP DOWN lag predicted by ripples.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (150ms windows UP lag control)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (full windows UP lag control)\V1 delta power predicted by ripples.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (full windows UP lag control)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression (full windows)\V1 delta power predicted by ripples.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression (full windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all



hFig = openfig('P:\corticohippocampal_replay\V1-HPC bilateral interaction\mixed effect regression 250ms ripples (full windows)\Ripples predicted by DOWN UP delta power.fig', 'visible');
print(hFig, sprintf('P:/corticohippocampal_replay/V1-HPC bilateral interaction/mixed effect regression 250ms ripples (full windows)/%s.svg',hFig.Name), '-dsvg', '-painters');
close all

hFig = gcf;
print(hFig, sprintf('P:/PhD/Dissertation/Figure materials/M24064 20241212/KDE_ripple_V1_reactivationSleepChronic/Probe2/%s.svg',hFig.Name), '-dsvg', '-painters');
close all


hFig = gcf;
print(hFig, sprintf('P:/PhD/Dissertation/Figure materials/M24017_20240604/Probe2/%s.svg',hFig.Name), '-dsvg', '-painters');
close all



save_all_figures('P:\corticohippocampal_replay\V1-HPC sleep reactivation\context-selective ripple PSTH',[])

save_all_figures('P:\corticohippocampal_replay\V1-HPC behaviour',[])

save_all_figures('P:\corticohippocampal_replay\V1-HPC behaviour\M24016 Ephys Day 5',[])

