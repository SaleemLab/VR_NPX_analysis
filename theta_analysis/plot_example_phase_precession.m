function plot_example_phase_precession

folders = {'X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Q-BLU\Q-BLU_Day8_16x8','X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\N-BLU\N-BLU_Day5_16x4',...
    'X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Q-BLU\Q-BLU_Day6_16x3','X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Q-BLU\Q-BLU_Day4_16x2',...
    'X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Q-BLU\Q-BLU_Day7_16x1'};

all_indx = [71,85,77,70,20,10];%last is T1
dires = [2,2,1,1,1,2]; %last is T1
%indx_chosen
figure
PP = plotting_parameters;
count  = 1;

for j = 1 : length(folders)
    
    cd(folders{j})
    load('extracted_phase_precession_absolute_location.mat')
    load('extracted_position.mat')
   
    if j == 5
        idcs = all_indx([j j+1]);
        dirs = dires([j j+1]);
    else
        idcs = all_indx(j);
        dir = dires(j);
    end
    
    if length(idcs) > 1
        tracks = [2,1];
    else
        track=2;
    end
        
    for i = 1 : length(idcs)
        if length(idcs) > 1 & i == 1
            track = tracks(i);
            dir = dirs(i);
        elseif length(idcs) > 1 & i == 2
            track = tracks(i);
            dir = dirs(i);
        end
        track_linear = position.linear(track).linear(~isnan(position.linear(track).linear));
        track_times = position.linear(track).timestamps;
        track_x = position.x(position.linear(track).clean_track_Indices);
        track_y = position.y(position.linear(track).clean_track_Indices);
        
        subplot(3,4,count)
        plot([TPP(track).(sprintf('%s','spike_positions_',num2str(dir))){idcs(i)},TPP(track).(sprintf('%s','spike_positions_',num2str(dir))){idcs(i)}],...
            [rad2deg(TPP(track).(sprintf('%s','spike_phases_',num2str(dir))){idcs(i)})+180,rad2deg(TPP(track).(sprintf('%s','spike_phases_',num2str(dir))){idcs(i)})+180+360],...
            'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3,'LineStyle','none')
        xlabel('Linearized position')
        ylabel('Phase (degrees)')
        set(findall(gcf,'-property','FontSize'),'FontSize',16)
        
        % create colour gradient for 360 degrees
        n=360;
        red = (1/ceil(n/2)):(1/ceil(n/2)):1;
        nSones = ones(ceil(n/2),1);
        red = [red' ; nSones];
        blue = flipud(red);
        green = zeros(1,size(red,1));
        cols = [red green' blue];
        
        %% plot spikes over animal's trajectory, with colour representing phase
        subplot(3,4,count+1)
        plot(track_x,track_y,'o','MarkerSize',2,...
            'MarkerEdgeColor','none','MarkerFaceColor',[0.3 0.3 0.3])
        spikes_x = TPP(track).(sprintf('%s','spike_x_',num2str(dir))){idcs(i)};
        spikes_y = TPP(track).(sprintf('%s','spike_y_',num2str(dir))){idcs(i)};
        spikes_pase = TPP(track).(sprintf('%s','spike_phases_',num2str(dir))){idcs(i)};
        hold on
        for n=1:length(spikes_pase)
            hold on
            plot(spikes_x(n),spikes_y(n),'o','MarkerSize',5,...
                'MarkerEdgeColor','none','MarkerFaceColor',cols(ceil(rad2deg(spikes_pase(n)))+180,:))
        end
        c = colorbar;
        cmap = colormap(cols);
        c.Limits = [0 1];
        c.Label.String  = 'Phase';
        c.Ticks = [0,1];
        c.TickLabels = {'0','360'};
        set(findall(gcf,'-property','FontSize'),'FontSize',16)
        
        count = count+2;
    end
    
    
end


end