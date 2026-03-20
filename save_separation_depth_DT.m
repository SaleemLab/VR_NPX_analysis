mouse = 'M25026';
date = '20250208';
separation_depth = struct();

% for MEC - HVA
separation_depth.MEC_HVA = 3000;
separation_depth.V1_HPC = 1625;
separation_depth.HPC_thalamus = 500;
% best_channels = cell(1,2);
% best_channels{1,1}.surface_depth = [NaN,4000,4000,NaN,4000,NaN,4000,NaN];
% best_channels{1,1}.HVA_depth = [NaN,3500,3500,NaN,3500,NaN,3500,NaN];
% best_channels{1,1}.MEC_entry_depth = [NaN,2700,2700,NaN,2700,NaN,2700,NaN];
% best_channels{1,1}.MEC_theta_depth = [NaN,2000,2000,NaN,2000,NaN,2000,NaN];
% best_channels{1,1}.MEC_ripple_depth = [NaN,1500,1500,NaN,1500,NaN,1500,NaN];
% best_channels{1,1}.xcoord = [27,59,277,309,527,559,777,809];
% 
% best_channels{1,2}.surface_depth = [NaN,3700,3800,NaN,3900,NaN,4000,NaN];
% best_channels{1,2}.L4_depth = [NaN,6180,6060,NaN,5820,NaN,5895,NaN];
% best_channels{1,2}.L5_depth = [NaN,6180,6060,NaN,5820,NaN,5895,NaN];
% best_channels{1,2}.CA1_depth = [NaN,6180,6060,NaN,5820,NaN,5895,NaN];
% best_channels{1,2}.xcoord = [27,59,277,309,527,559,777,809];

save(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',mouse,'\analysis\','separation_depth.mat'),"separation_depth");