%% Tomaso Muzzu - UCL - 12/12/2018

%% version for Bonvision stimulus

%% photodiode signal: save original signal downsampled at 1kHz, timestamps of stimulus onsets and offsets
function [ACInfo] = OE_Analyse_PD_Signal(ACInfo,frame_dur_tol)
    
%     if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
%         load([DIRname filesep 'AC_Info.mat'],'stimON','stimOFF');
%         % else load the analog channels now, find when the timestamps of the VS
%         % onsets and offsets
%     else
%         display('Analogue channels have not been converted yet or some error reading the saved file. Please select folder containing relevant OE recording for recovering ACs info.');
        
        % data(:,1) = sync pulse signal
        % data(:,2) = photodiode signal
        % data(:,3) = signal A from rotary encoder
        % data(:,4) = signal B from rotary encoder
        % data(:,5) = Bonvision state signal      
        
        % resample and normalise signal
        if size(ACInfo.Data,2)>5
            normPD = ACInfo.Data(:,4)/max(ACInfo.Data(:,4));
        else
            normPD = ACInfo.Data(:,2)/max(ACInfo.Data(:,2));
        end
        [temp_pks, temp_locs] = findpeaks(normPD);
        GMModel = fitgmdist(normPD,2);
        ThresholdPeaks = max(GMModel.mu);
        pks = temp_pks(temp_pks>ThresholdPeaks);
        locs = temp_locs(temp_pks>ThresholdPeaks);
        clear temp_pks temp_locs
%         figure
%         plot(plot_TS_temp,plot_PD_temp)
%         hold on
%         plot(ACInfo.Timestamps(locs),pks,'.')
        
        
        % find first peaks signalling the onset of the white square
        PksTimeDifference = diff(locs);
        clear temp_VSOnsets VSOnsetsIndecesOFF VSOnsetsIndecesON
        temp_VSOnsets = find(PksTimeDifference>frame_dur_tol*(ACInfo.SamplingRateOE) & PksTimeDifference<0.15*(ACInfo.SamplingRateOE)); 
        VSOnsetsIndeces = sort([locs(temp_VSOnsets); locs(temp_VSOnsets+1)]); % save indexes of stimulus onsets

        % first change of the sync square
        downsample_ratio = 15;
        plot_TS_temp = resample(ACInfo.Timestamps,1,downsample_ratio);
        plot_PD_temp = resample(normPD,1,downsample_ratio);
        
        recDate = str2num([ACInfo.ExpDate(1:4) ACInfo.ExpDate(6:7) ACInfo.ExpDate(9:10)]);
        if isempty(recDate)
            recDate = str2num([ACInfo.ExpDate(1:2) ACInfo.ExpDate(4:5) ACInfo.ExpDate(7:8)]);
            recDate = 20000000+recDate;
        end
        if recDate < 20190819 & recDate > 20190101
            % look at the PD signal preceding the first PD change found above
            % initial gray screen has white sync square
            if sum(normPD(VSOnsetsIndeces(1)-1*ACInfo.SamplingRateOE:VSOnsetsIndeces(1))>ThresholdPeaks*0.6) >= ...
                    length(normPD(VSOnsetsIndeces(1)-1*ACInfo.SamplingRateOE:VSOnsetsIndeces(1)))/4
                ACInfo.SE_periods(1,1:2) = [ACInfo.Timestamps(VSOnsetsIndeces(1))-mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))-1,...
                    ACInfo.Timestamps(VSOnsetsIndeces(1))-mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))];
                [v_ ind_1] = min(abs(ACInfo.Timestamps-ACInfo.SE_periods(1,2)));
                VSOnsetsIndeces = [ind_1 ; VSOnsetsIndeces];
            else % initial gray screen has black sync square
                ACInfo.SE_periods(1,1:2) = [ACInfo.Timestamps(VSOnsetsIndeces(1))-2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))-1,...
                                            ACInfo.Timestamps(VSOnsetsIndeces(1))-2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))];
                [v_ ind_1] = min(abs(ACInfo.Timestamps-(ACInfo.SE_periods(1,2)+0.1)));
                [v_ ind_2] = min(abs(ACInfo.Timestamps-ACInfo.SE_periods(1,2)));
                VSOnsetsIndeces = [ind_1; ind_2; VSOnsetsIndeces];
            end
            % final gray screen has white sync square
            if  sum(normPD(VSOnsetsIndeces(end):VSOnsetsIndeces(end)+1*ACInfo.SamplingRateOE)>ThresholdPeaks*0.6) >= ...
                    length(normPD(VSOnsetsIndeces(end):VSOnsetsIndeces(end)+1*ACInfo.SamplingRateOE))/4
                ACInfo.SE_periods(2,1:2) = [ACInfo.Timestamps(VSOnsetsIndeces(end))+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))),...
                                            ACInfo.Timestamps(VSOnsetsIndeces(end))+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))+1];
                [v_ ind_1] = min(abs(ACInfo.Timestamps-ACInfo.SE_periods(2,1)));
                VSOnsetsIndeces = [ind_1; VSOnsetsIndeces];
            else % final gray screen has black sync square
                ACInfo.SE_periods(2,1:2) = [ACInfo.Timestamps(VSOnsetsIndeces(end))+2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))),...
                                            ACInfo.Timestamps(VSOnsetsIndeces(end))+2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))+1];
                [v_ ind_1] = min(abs(ACInfo.Timestamps-(ACInfo.SE_periods(2,1)-0.1)));
                [v_ ind_2] = min(abs(ACInfo.Timestamps-ACInfo.SE_periods(2,1)));
                VSOnsetsIndeces = [ind_1; ind_2 ; VSOnsetsIndeces];                       
            end
            ACInfo.PD_changes = unique(sort(ACInfo.Timestamps(VSOnsetsIndeces)));
        elseif recDate < 20190101
            ACInfo.PD_changes = unique(sort([ACInfo.Timestamps(VSOnsetsIndeces); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(1))-mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(end))+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(end))+2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))]));
            ACInfo.SE_periods(1,1:2) = [min(ACInfo.PD_changes)-1 min(ACInfo.PD_changes)];
            ACInfo.SE_periods(2,1:2) = [max(ACInfo.PD_changes) max(ACInfo.PD_changes)+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))+1];
        elseif recDate > 20190819
            ACInfo.PD_changes = unique(sort([ACInfo.Timestamps(VSOnsetsIndeces); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(1))-mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(end))+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces))); ...
                                             ACInfo.Timestamps(VSOnsetsIndeces(end))+2*mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))]));
            ACInfo.SE_periods(1,1:2) = [min(ACInfo.PD_changes)-1 min(ACInfo.PD_changes)];
            ACInfo.SE_periods(2,1:2) = [max(ACInfo.PD_changes) max(ACInfo.PD_changes)+mean(diff(ACInfo.Timestamps(VSOnsetsIndeces)))+1];
        end
        
        % verify the events are detected correctly:
        fprintf('Please verify that all trial and perturbation onsets/offsets are detected correctly\n');
        figure
        plot(plot_TS_temp,plot_PD_temp)
        hold on
        plot(ACInfo.PD_changes,0.7*ones(1,length(ACInfo.PD_changes)),'.r','MarkerSize',18)
        hold on % initial period
        plot(ACInfo.SE_periods(1,1),0.5*ones(1,1),'ok','MarkerSize',5)
        plot(ACInfo.SE_periods(1,2),0.5*ones(1,1),'or','MarkerSize',8)
        hold on % final period
        plot(ACInfo.SE_periods(2,1),0.5*ones(1,1),'ok','MarkerSize',5)
        plot(ACInfo.SE_periods(2,2),0.5*ones(1,1),'or','MarkerSize',8)

end    
