function plot_ripple_power_UP(slow_waves_all,ripples_all,event_info,probability_psth,probability_normalised,probability_psth_whole,probability_normalised_whole,UP_DOWN_ripple_PSTH_MUA)

if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

%% Select UP events based on ripple power
bin_edges = [0 0.2 0.4 0.6 0.8 1];

for nprobe = 1:2
   
    % 0 0.2 0.4 0.6 0.8
    for mprobe = 1:2
         SWS_ripples_index = find(ripples_all(mprobe).SWS_index);
        for nbin = 1:length(bin_edges)-1
            UP_index{nprobe}{mprobe}{nbin} = []; % index out of selected UP events for grabbing already-calculated metrics
            UP_index_all{nprobe}{mprobe}{nbin} = []; % index out of all UP events (can be useful for grabbing additional information not previously calculated)
            DOWN_index{nprobe}{mprobe}{nbin} = []; % index out of selected DOWN events for grabbing already-calculated metrics
            DOWN_index_all{nprobe}{mprobe}{nbin} = []; % index out of all DOWN events (can be useful for grabbing additional information not previously calculated)

            ripples_index{nprobe}{mprobe}{nbin} = [];% index out of all ripple events (can be useful for grabbing additional information not previously calculated) 
            ripples_index_all{nprobe}{mprobe}{nbin} = [];% index out of selected ripple events (can be useful for grabbing additional information not previously calculated) 

            if mprobe == 1
                UP_index_all{nprobe}{mprobe}{nbin} = event_info(nprobe).L_ripple_normalised_UP_duration(event_info(nprobe).L_ripple_normalised_UP_duration(:,3) > bin_edges(nbin) & event_info(nprobe).L_ripple_normalised_UP_duration(:,3) < bin_edges(nbin+1),2);
                [C,ia,ib] = intersect(event_info(nprobe).UP_index,UP_index_all{nprobe}{mprobe}{nbin});
                UP_index{nprobe}{mprobe}{nbin} = ia;

                DOWN_index_all{nprobe}{mprobe}{nbin} = event_info(nprobe).L_ripple_normalised_DOWN_duration(event_info(nprobe).L_ripple_normalised_DOWN_duration(:,3) > bin_edges(nbin) & event_info(nprobe).L_ripple_normalised_DOWN_duration(:,3) < bin_edges(nbin+1),2);
                [C,ia,ib] = intersect(event_info(nprobe).DOWN_index,DOWN_index_all{nprobe}{mprobe}{nbin});
                DOWN_index{nprobe}{mprobe}{nbin} = ia;

                ripples_index_all{nprobe}{mprobe}{nbin} = event_info(nprobe).L_ripple_normalised_UP_duration(event_info(nprobe).L_ripple_normalised_UP_duration(:,3) > bin_edges(nbin) & event_info(nprobe).L_ripple_normalised_UP_duration(:,3) < bin_edges(nbin+1),1);
            else
                UP_index_all{nprobe}{mprobe}{nbin} =  event_info(nprobe).R_ripple_normalised_UP_duration(event_info(nprobe).R_ripple_normalised_UP_duration(:,3) > bin_edges(nbin) & event_info(nprobe).R_ripple_normalised_UP_duration(:,3) < bin_edges(nbin+1),2);
                [C,ia,ib] = intersect(event_info(nprobe).UP_index,UP_index_all{nprobe}{mprobe}{nbin});
                UP_index{nprobe}{mprobe}{nbin} = ia;

                DOWN_index_all{nprobe}{mprobe}{nbin} = event_info(nprobe).R_ripple_normalised_DOWN_duration(event_info(nprobe).R_ripple_normalised_DOWN_duration(:,3) > bin_edges(nbin) & event_info(nprobe).R_ripple_normalised_DOWN_duration(:,3) < bin_edges(nbin+1),2);
                [C,ia,ib] = intersect(event_info(nprobe).DOWN_index,DOWN_index_all{nprobe}{mprobe}{nbin});
                DOWN_index{nprobe}{mprobe}{nbin} = ia;

                ripples_index_all{nprobe}{mprobe}{nbin} = event_info(nprobe).R_ripple_normalised_UP_duration(event_info(nprobe).R_ripple_normalised_UP_duration(:,3) > bin_edges(nbin) & event_info(nprobe).R_ripple_normalised_UP_duration(:,3) < bin_edges(nbin+1),1);
            end

            [C,ia,ib] = intersect(SWS_ripples_index,ripples_index_all{nprobe}{mprobe}{nbin});
            ripples_index{nprobe}{mprobe}{nbin} = ia;
        end
    end
end

%%

colour_lines = [215,48,39;244,109,67;145,191,219;69,117,180]/256;% 4 colors Dark red, orangish red, light blue, dark blue
clear colour_lines

colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue 

%% Plotting distribution of ripple during normalised duration of UP  (peaktimes)
probability = probability_normalised;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,48,39;244,109,67;145,191,219;69,117,180]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue 
probe_hemisphere_texts = {'Probability of ripples at different timing during left V1 normalised UP-DOWN duration','Probability of ripples at different timing during right V1 normalised UP-DOWN duration'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
probability_ripple_timing_normalised = [];
% colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);
    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).L_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_DOWN(index,:),1),size(probability(nprobe).L_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).L_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).R_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_DOWN(index,:),1),size(probability(nprobe).R_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).R_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% Left ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).L_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_UP(index,:),1),size(probability(nprobe).L_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).L_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).R_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_UP(index,:),1),size(probability(nprobe).R_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).R_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end



%% plotting Probability of SWR relative to UP and DOWN onset (peaktimes)
probability = probability_psth;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,48,39;244,109,67;145,191,219;69,117,180]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue
probe_hemisphere_texts = {'Probability of ripples at different timing during left V1 UP-DOWN','Probability of ripples at different timing during right V1 UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
probability_ripple_timing_normalised = [];
% colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);
    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).L_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_DOWN(index,:),1),size(probability(nprobe).L_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).L_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).R_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_DOWN(index,:),1),size(probability(nprobe).R_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).R_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% Left ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).L_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_UP(index,:),1),size(probability(nprobe).L_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).L_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).R_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_UP(index,:),1),size(probability(nprobe).R_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).R_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end


if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save(fullfile(analysis_folder,'V1-HPC sleep interaction','probability_ripple_timing_normalised.mat'),'probability_ripple_timing_normalised');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','probability_ripple_timing.mat'),'probability_ripple_timing');

%% Plotting distribution of ripple during normalised duration of UP  (WHOLE)
probability = probability_normalised_whole;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,48,39;244,109,67;145,191,219;69,117,180]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue 
probe_hemisphere_texts = {'Probability of WHOLE ripples at different timing during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples at different timing during right V1 normalised UP-DOWN duration'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
probability_ripple_timing_normalised = [];
% colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);
    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).L_ripples_DOWN(index,:),2))/all_DOWN_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_DOWN(index,:),1),size(probability(nprobe).L_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_DOWN(index(event_id),:))/all_DOWN_no;
        end
        probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).L_ripples_DOWN(index,:))/all_DOWN_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).R_ripples_DOWN(index,:),2))/all_DOWN_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_DOWN(index,:),1),size(probability(nprobe).R_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_DOWN(index(event_id),:))/all_DOWN_no;
        end
        probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of DOWN')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).R_ripples_DOWN(index,:))/all_DOWN_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Normalised duration of DOWN')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% Left ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).L_ripples_UP(index,:),2))/all_UP_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_UP(index,:),1),size(probability(nprobe).L_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_UP(index(event_id),:))/all_UP_no;
        end
        probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).L_ripples_UP(index,:))/all_UP_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).L_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(cumsum(probability(nprobe).R_ripples_UP(index,:),2))/all_UP_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_UP(index,:),1),size(probability(nprobe).R_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_UP(index(event_id),:))/all_UP_no;
        end
        probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);
        
        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = linspace(0,1,num_bins);
        y = sum(probability(nprobe).R_ripples_UP(index,:))/all_UP_no;
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing_normalised(nprobe).R_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Normalised duration of UP')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end



%% plotting Probability of SWR relative to UP and DOWN onset (WHOLE)
probability = probability_psth_whole;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,48,39;244,109,67;145,191,219;69,117,180]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue
probe_hemisphere_texts = {'Probability of WHOLE ripples at different timing during left V1 UP-DOWN','Probability of WHOLE ripples at different timing during right V1 UP-DOWN'};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 normalised UP-DOWN duration','Probability of WHOLE ripples during right V1 normalised UP-DOWN duration'};
probability_ripple_timing_normalised = [];
% colour_lines = [215,25,28;44,123,182]/256;
for nprobe = 1:2
    fig(nprobe)=figure;
    fig(nprobe).Position = [982 50 700 950];
    fig(nprobe).Name = probe_hemisphere_texts{nprobe};

    % all_ripple_no = probability(nprobe).L_ripple_no;
    all_UP_no = length(probability(nprobe).UP_all_index);
    all_DOWN_no = length(probability(nprobe).DOWN_all_index);
    %%%% Left ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).L_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_DOWN(index,:),1),size(probability(nprobe).L_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).L_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).L_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% DOWN
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).R_ripples_DOWN(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_DOWN(index,:),1),size(probability(nprobe).R_ripples_DOWN(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_DOWN(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).R_ripples_DOWN(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).R_ripple_DOWN_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during DOWN')
    xlabel('Time relative to DOWN onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


    %%%% Left ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).L_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).L_ripples_UP(index,:),1),size(probability(nprobe).L_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).L_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).L_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).L_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



    %%%% Right ripples plotting
    %%%% UP
    nexttile
    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(cumsum(probability(nprobe).R_ripples_UP(index,:),2))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        tempUP = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
            event_id = datasample(s,1:size(probability(nprobe).R_ripples_UP(index,:),1),size(probability(nprobe).R_ripples_UP(index,:),1));
            tempUP(iBoot,:) = sum(probability(nprobe).R_ripples_UP(index(event_id),:))/length(index(event_id));
        end
        probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin} = tempUP;

        LCI = prctile(cumsum(tempUP,2),2.5);
        UCI = prctile(cumsum(tempUP,2),97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    % title('Probability of left ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Cumulative probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    nexttile

    for nbin = 1:5
        index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));

        x = time_wondows(1)+time_bin/2:time_bin:time_wondows(end)-time_bin/2;
        y = sum(probability(nprobe).R_ripples_UP(index,:))/length(index);
        %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

        LCI = prctile(probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin},2.5);
        UCI = prctile(probability_ripple_timing(nprobe).R_ripple_UP_bootstrap{nbin} ,97.5);

        PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
        ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
    end
    % legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
    % xline(0,'r')
    title('Probability of right ripples during UP')
    xlabel('Time relative to UP onset (s)')
    ylabel('Probability')
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
save(fullfile(analysis_folder,'V1-HPC sleep interaction','probability_ripple_timing_normalised_whole.mat'),'probability_ripple_timing_normalised');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','probability_ripple_timing_whole.mat'),'probability_ripple_timing');


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])




%% MUA during UP with different ripple timing
% probability = probability_psth;
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue

% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left V1 UP with different left ripple timing','MUA activity during right V1 UP with different left ripple timing',...
    'MUA activity during left V1 UP with different right ripple timing','MUA activity during right V1 UP with different right ripple timing',};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
ncount = 1;
clear MUA_ripple_timing
for mprobe = 1:2
    for nprobe = 1:2
        fig(ncount)=figure;
        fig(ncount).Position = [1355 332 700 600];
        fig(ncount).Name = probe_hemisphere_texts{ncount};
        ncount = ncount+ 1;

        all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
        % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

        %%%% Left V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).L_V1_UP(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).L_V1_UP(index,:),1),size(PSTH_MUA(nprobe).L_V1_UP(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).L_V1_UP(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).L_V1_UP_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Left V1')
        xlabel('Time relative to UP onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        %%% Right V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).R_V1_UP(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).R_V1_UP(index,:),1),size(PSTH_MUA(nprobe).R_V1_UP(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).R_V1_UP(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_V1_UP_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Right V1')
        xlabel('Time relative to UP onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



        %%%% Left HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).L_HPC_UP(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).L_HPC_UP(index,:),1),size(PSTH_MUA(nprobe).L_HPC_UP(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).L_HPC_UP(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).L_HPC_UP_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
        % xline(0,'r')
        title('Left HPC')
        xlabel('Time relative to UP onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        %%% Right HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).R_HPC_UP(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).R_HPC_UP(index,:),1),size(PSTH_MUA(nprobe).R_HPC_UP(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).R_HPC_UP(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_HPC_UP_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
        % xline(0,'r')
        title('Right HPC')
        xlabel('Time relative to UP onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
end






%% MUA during ripples with different ripple timing
% probability = probability_psth;
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue

% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left ripple with different timing relative to left V1','MUA activity during left ripple with different timing relative to right V1',...
    'MUA activity during right ripple with different timing relative to left V1','MUA activity during right ripple with different timing relative to right V1',};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
ncount = 1;
for mprobe = 1:2
    for nprobe = 1:2
        fig(ncount)=figure;
        fig(ncount).Position = [1355 332 700 600];
        fig(ncount).Name = probe_hemisphere_texts{ncount};
        ncount = ncount+ 1;

        all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
        % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

        %%%% Left V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,ripples_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(mprobe).L_V1_ripples(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(mprobe).L_V1_ripples(index,:),1),size(PSTH_MUA(mprobe).L_V1_ripples(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(mprobe).L_V1_ripples(index(event_id),:));
            end
            MUA_ripple_timing(mprobe).L_V1_ripples_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Left V1')
        xlabel('Time relative to ripple peaktimes (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        %%% Right V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,ripples_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(mprobe).R_V1_ripples(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(mprobe).R_V1_ripples(index,:),1),size(PSTH_MUA(mprobe).R_V1_ripples(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(mprobe).R_V1_ripples(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_V1_ripples_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Right V1')
        xlabel('Time relative to ripple peaktimes (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



        %%%% Left HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,ripples_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(mprobe).L_HPC_ripples(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(mprobe).L_HPC_ripples(index,:),1),size(PSTH_MUA(mprobe).L_HPC_ripples(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(mprobe).L_HPC_ripples(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).L_HPC_ripples_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
        % xline(0,'r')
        title('Left HPC')
        xlabel('Time relative to ripple peaktimes (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        %%% Right HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,ripples_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(mprobe).R_HPC_ripples(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(mprobe).R_HPC_ripples(index,:),1),size(PSTH_MUA(mprobe).R_HPC_ripples(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(mprobe).R_HPC_ripples(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_HPC_ripples_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
        % xline(0,'r')
        title('Right HPC')
        xlabel('Time relative to ripple peaktimes (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
end


%% MUA during DOWN with different ripple timing
% probability = probability_psth;
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;

time_wondows = [-0.2 0.5];
time_bin = 0.02;
num_bins=20; % divide one UP event into 20 bins
duration_threshold = 2;
all_sessions = max(slow_waves_all(1).DOWN_session_count);
% colour_lines = [215,25,28;214,100,77;103,169,207;44,123,182]/256;
clear colour_lines
colour_lines{1} = [252,146,114;251,106,74;239,59,44;203,24,29;153,0,13]/256;% 5 red
colour_lines{2} = [158,202,225;107,174,214;66,146,198;33,113,181;8,69,148]/256;% 5 blue

% colour_lines = [215,25,28;44,123,182]/256;
probe_hemisphere_texts = {'MUA activity during left V1 DOWN with different left ripple timing','MUA activity during right V1 DOWN with different ripple timing',...
    'MUA activity during left V1 DOWN with different right ripple timing','MUA activity during right V1 DOWN with different right ripple timing',};
% probe_hemisphere_texts = {'Probability of WHOLE ripples during left V1 UP-DOWN','Probability of WHOLE ripples during right V1 UP-DOWN'};
ncount = 1;
clear MUA_ripple_timing
for mprobe = 1:2
    for nprobe = 1:2
        fig(ncount)=figure;
        fig(ncount).Position = [1355 332 700 600];
        fig(ncount).Name = probe_hemisphere_texts{ncount};
        ncount = ncount+ 1;

        all_UP_no = size(PSTH_MUA(nprobe).L_V1_DOWN,1);
        % all_DOWN_no = length(probability(nprobe).DOWN_all_index);

        %%%% Left V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,UP_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).L_V1_DOWN(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).L_V1_DOWN(index,:),1),size(PSTH_MUA(nprobe).L_V1_DOWN(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).L_V1_DOWN(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).L_V1_DOWN_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Left V1')
        xlabel('Time relative to DOWN onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

        %%% Right V1
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,DOWN_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).R_V1_DOWN(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).R_V1_DOWN(index,:),1),size(PSTH_MUA(nprobe).R_V1_DOWN(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).R_V1_DOWN(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_V1_DOWN_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end

        title('Right V1')
        xlabel('Time relative to DOWN onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



        %%%% Left HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,DOWN_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).L_HPC_DOWN(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).L_HPC_DOWN(index,:),1),size(PSTH_MUA(nprobe).L_HPC_DOWN(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).L_HPC_DOWN(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).L_HPC_DOWN_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        % legend(ERROR_SHADE(1:2),{'Original','Shuffled'},'Box','off')
        % xline(0,'r')
        title('Left HPC')
        xlabel('Time relative to DOWN onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


        %%% Right HPC
        nexttile
        for nbin = 1:5
            tic
            index = unique(cat(1,DOWN_index{nprobe}{mprobe}{nbin}));
            x = PSTH_MUA(nprobe).timebins;
            y = mean(PSTH_MUA(nprobe).R_HPC_DOWN(index,:));
            %     y = mean(cumsum(probability(nprobe).L_ripples_DOWN_bootstrap,2));

            tempUP = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
                event_id = datasample(s,1:size(PSTH_MUA(nprobe).R_HPC_DOWN(index,:),1),size(PSTH_MUA(nprobe).R_HPC_DOWN(index,:),1));
                tempUP(iBoot,:) = mean(PSTH_MUA(nprobe).R_HPC_DOWN(index(event_id),:));
            end
            MUA_ripple_timing(nprobe).R_HPC_DOWN_bootstrap{mprobe}{nbin} = tempUP;

            LCI = prctile(tempUP,2.5);
            UCI = prctile(tempUP,97.5);
            % LCI = y+ (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));
            % UCI = y- (std(PSTH_MUA(nprobe).L_V1_UP(index,:))/length(index));

            PLOT = plot(x,y,'Color',colour_lines{nprobe}(nbin,:));hold on;
            ERROR_SHADE(nbin) = patch([x fliplr(x)],[UCI fliplr(LCI)],colour_lines{nprobe}(nbin,:),'FaceAlpha','0.3','LineStyle','none');
            xlim([-0.5 0.5])
            toc
        end
        legend(ERROR_SHADE(1:5),{'0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'},'Box','off')
        % xline(0,'r')
        title('Right HPC')
        xlabel('Time relative to DOWN onset (s)')
        ylabel('MUA activty (z)')
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
end

save_all_figures(fullfile(analysis_folder,'V1-HPC sleep interaction'),[])

if exist(fullfile(analysis_folder,'V1-HPC sleep interaction')) ==0
    mkdir(fullfile(analysis_folder,'V1-HPC sleep interaction'))
end
% save(fullfile(analysis_folder,'V1-HPC sleep interaction','probability_ripple_timing_normalised_whole.mat'),'probability_ripple_timing_normalised');
save(fullfile(analysis_folder,'V1-HPC sleep interaction','MUA_ripple_timing.mat'),'MUA_ripple_timing');



