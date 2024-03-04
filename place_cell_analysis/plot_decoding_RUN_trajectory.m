function plot_decoding_RUN_trajectory(estimated_position_lap_CV,Task_info,options)


pcount = 1;
nfigure = 1;
for track_id = 1:2
    track_laps = find(Task_info.track_ID_all == track_id);
    [C,ia,ib] = intersect(Task_info.complete_laps_id,track_laps);

    track_laps = ib;

    for lap_id =round(linspace(1,length(track_laps),16))


        %             for lap_id = lap_times(track_id).completeLaps_id(6:2:40)
        if pcount == 17
            nfigure = nfigure + 1;
            pcount = 1;
        end

        fig = figure(nfigure)
        fig.Position = [300 150 945 800];
        fig.Name = (sprintf('%s %s CV %s Bayesian decoding visualisation',options.SUBJECT,options.SESSION,options.region));
        subplot(4,4,pcount)
        if ~isempty(estimated_position_lap_CV(track_id).lap(lap_id))
            imagesc([estimated_position_lap_CV(track_id).lap(lap_id).track(1).run; estimated_position_lap_CV(track_id).lap(lap_id).track(2).run])
            hold on
            plot(estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_actual_position/10,'r')
            plot(estimated_position_lap_CV(track_id).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
            yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
            yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
            yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
            run_time_edges = estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_time_edges;

            xticks(linspace(1,length(run_time_edges),5))
            xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')
        title(sprintf('Lap %i',lap_id))
        colorbar
        colormap(flip(bone))
        pcount = pcount + 1;
    end
end

sgtitle((sprintf('%s %s CV %s Bayesian decoding visualisation',options.SUBJECT,options.SESSION,options.region)))