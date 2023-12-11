classdef NPfra
    %UNTITLED Functions to plot FRAs for NP
    %   Detailed explanation goes here
    
    methods (Static)
        function plotSingleRec(binName,binPath)
            % check if its been sorted by Phy
            sortedCheck = dir([binPath,'kilosort3\','phy.log']);
            if isempty(sortedCheck)
                errordlg('Not sorted yet')
                return
            end
            
            % FRA params
            binWidth = 0.01; % In seconds
            ISI = 0.1;
            StimLength = 0.1;
            bHvDir = 'C:\Users\Kat\Dropbox\Data\';
            % Get ferret name
            filepathparts = strsplit(binPath,'\');
            currFerret = filepathparts{min(find(contains(filepathparts,'F')))};
            % Get NP trial times
            [trialTimes,~] = NPadmin.getTrialTimes(binName,binPath);
            % Get meta
            [meta] = NPadmin.ReadMeta(binName, binPath);
            
            % Get behavioural FRA data
            currNPDateNum = addtodate(datenum(meta.fileCreateTime,'yyyy-mm-ddTHH:MM:SS'),round(meta.fileTimeSecs),'second');
            bHvFiles = dir([bHvDir,currFerret,'\*log.mat']);
            % this is the date modified so isnt' created
            allDateNum = cell2mat({bHvFiles.datenum});
            [miny,mini] = min(abs(allDateNum - currNPDateNum)); % Finds smallest time
            bHvfileName = bHvFiles(mini).name;
            load([bHvDir,currFerret,'\',bHvfileName])
            
            if numel(trialTimes)~=numel(StartSamples)
                [trialTimes2,StartSamples2] = alignsignals(trialTimes,StartSamples);
                % then remove zeroes 
                StartSamples3 = StartSamples2(StartSamples2~=0);
                trialTimes3 = trialTimes2(StartSamples2~=0);
                trialTimes = trialTimes3(trialTimes3 ~=0);
                StartSamples = StartSamples3(trialTimes3 ~=0);
                
                disp('Bhv and NI trial times dont equal each other. Check plot')
                figure;plot(diff(StartSamples),diff(trialTimes));
            end
            % Get spike data
            spike = NPksAdmin.extractPhySpikes([binPath,'kilosort3\']);
            
            % Associate spike times with spike Sites
            spikeData(:,1) = spike.Times;
            spikeData(:,2) = spike.Sites;
            % Associate with depth
            spikeData(:,3) = spike.Ycoords(spike.Sites);
            % Associate with shank
            
            Shanks = NPadmin.xcoordsToShank(spike.Xcoords);
            spikeData(:,4) = Shanks(spike.Sites);
            % Get metadata
            freqs = unique(cell2mat(passiveRegStim(:,2)));
            nFreqs = length(freqs);
            levels = unique(cell2mat(passiveRegStim(:,3)));
            nLevel = length(levels);
            sitedepths = sortrows(unique(spikeData(:,2:4),'rows'),[3 2 1]);
            
            nChan = size(sitedepths,1);
            
            spikeRateSpont = NaN(length(freqs)*length(levels),nChan);
            spikeRateResp = NaN(length(freqs)*length(levels),nChan);
            
            % Treating a channel as a cluster at a unique site/depth
            for currSite = 1:nChan
                currSpikeTimes = spikeData(spikeData(:,3)==sitedepths(currSite,2)&spikeData(:,2)==sitedepths(currSite,1),1);
                
                % Calculate onset latency
                [~,~,~,ally]= NPplot.NPgetPSTH(currSpikeTimes,trialTimes,ISI,StimLength,binWidth);
                % Want a longer window for FRA PSTH
                [yfra,ts,~,allyfra]= NPplot.NPgetPSTH(currSpikeTimes,trialTimes,ISI,0.2,binWidth/2);
                PSTH(currSite,:) = yfra;
                PSTHse(currSite,:) = std(allyfra)/sqrt(size(allyfra,1));
                [clusters, p_values, t_sums, permutation_distribution ] = Stats.clusterStats(ally(:,(1:size(ally,2)/2))', ally(:,size(ally,2)/2+1:end)',0,0.05,1000,1,[]);
                disp([num2str(currSite),'/',num2str(nChan)]);
                
                % Remove clusters above the alpha
                clusters = clusters(p_values < 0.05);
                latencies = [];
                onlength = [];
                for i = 1:size(clusters,2)
                    latencies(i) = min(clusters{1,i})*binWidth; % in ms
                    onlength(i) = size(clusters{1,i},2);
                end
                if isempty(latencies)
                    OnsetLatencys(currSite) = NaN;
                    OnsetLengths(currSite) = NaN;
                else
                    OnsetLatencys(currSite) = min(latencies);
                    OnsetLengths(currSite) = sum(onlength)*binWidth;
                end
                % Get the mean PSTH for each trial idx
                for currTrialIdx = unique(trialIdx)
                    currTrialTimes = trialTimes(trialIdx==currTrialIdx);
                    [y,~]= NPplot.NPgetPSTH(currSpikeTimes,currTrialTimes,ISI,StimLength,binWidth);
                    % Just takes the halfway point
                    spikeRateSpont(currTrialIdx,currSite) = mean(y(1:length(y/2)));
                    spikeRateResp(currTrialIdx,currSite) = mean(y(length(y)/2+1:end));
                end
                
            end
            
            % Preassign
            [fra, n] = deal(nan(nLevel, nFreqs, nChan));
            
            % get numbers for FRA
            for i = 1 : nFreqs
                % For each level
                for j = 1 : nLevel
                    for currChan = 1:nChan
                        % Collate all trials
                        [~,Idx] = ismember([freqs(i) levels(j)],cell2mat(passiveRegStim(:,2:3)),'rows');
                        temp = spikeRateResp(Idx,currChan);
                        sponTemp = spikeRateSpont(Idx,currChan);
                        
                        % Average across trials
                        fra(j,i,currChan) = nanmean(temp);
                        n(j,i,currChan)   = numel(temp);
                        spon(j,i,currChan) = nanmean(sponTemp); % Need this for calculating CF
                    end
                end
            end
            
            % Convert freqs into kHz
            newfreqs = freqs ./ 1e3;
            freq_tick = 1 : 5 : nFreqs;
            freq_val  = num2str(round(newfreqs(freq_tick),2));
            
            yfiglim = 3; xfiglim = 10;
            totalfigs = ceil(nChan/(yfiglim*xfiglim));
            chan = 1;
            for f = 1:totalfigs
                fig(f) = figure('units','centimeters','position',[2 2 45 23],...
                    'name',sprintf('FRA_%d', f));
                sp = [];
                
                for i = 1 : yfiglim
                    for j = 1 : xfiglim
                        % Create axes
                        h = 1.5+((j-1)*4); % horizontal position
                        v = 3.9+((i-1)*7.1); % vertical position
                        
                        v2 = 1+((i-1)*7.1); % vertical position
                        sp(i,j) = axes('units','centimeters',...
                            'position',[h v 3 3]);
                        
                        z = fra(:,:,chan);
                        ss = spon(:,:,chan);
                        % Do some smoothing
                        z = NPfra.smoothFRA(z);
                        
                        % Since the largest two levels of the lowest
                        % frequences clip need to remove them
                        z(end,1:2) = [0 0]; % 0 for it to not be included in the bounds calculation
                        
                        spikes = z;
                        %calculate the mean spontaneous rate and calculate a threshold based on spontaneus
                        %rate+20% of max rate during stimulus-spontaneus rate.
                        
                        srate=mean(ss,'all')+(1/5*(max(max(spikes))));
                        
                        % calculate the bounds
                        [bounds,CF,Q10,Q30] = NPfra.calculateFRAbounds(spikes,nFreqs,nLevel,newfreqs.*1e3,srate);
                        %   Normalise
                        z(end,1:2) = [NaN NaN]; % Needs to be nan
                        
                        z = z - min(z(:));  % Subtract minimum
                        z = z ./ max(z(:)); % Normalize
                        
                        % Smoother and prettier than imagesc
                        contourf(1:nFreqs,levels,z,40,'linecolor','none'); hold on
                        
                        % Calculate CF
                        if ~isempty(bounds)
                            plot(1:nFreqs,bounds*10,'color','w','LineWidth',2); hold on
                            if ~isnan(CF)
                                % plot Q10
                                plot(find(ismember(newfreqs.*1e3,Q10)),[min(bounds+1)*10 min(bounds+1)*10],'r','LineWidth',1); hold on
                                % plot Q30
                                plot(find(ismember(newfreqs.*1e3,Q30)),[min(bounds+3)*10 min(bounds+3)*10],'m','LineWidth',1); hold on
                            end
                        end
                        CharactersticFrequency = round(CF);
                        
                        
                        title(['S',num2str(sitedepths(chan,end)),' Depth: ',num2str(sitedepths(chan,2))]);
                        if ~isnan(OnsetLatencys(chan))
                            xlabel(['* CF: ',num2str(CharactersticFrequency)]);
                        else
                            xlabel(['CF: ',num2str(CharactersticFrequency)]);
                        end
                        storeCF(chan) = CharactersticFrequency;
                        storeBounds(chan,:)=bounds;
                        storeQ10(chan,:) = Q10;
                        storeQ30(chan,:) = Q30;
                        ylim([levels(1) levels(end)])
                        
                        % Plot PSTH
                        pst(i,j) = axes('units','centimeters',...
                            'position',[h v2 3 1.5]);
                        patch([ts flip(ts)],[PSTH(chan,:)-PSTHse(chan,:) flip(PSTH(chan,:)+PSTHse(chan,:))],[.6 .6 .6],'EdgeAlpha',0,'FaceAlpha',.5); hold on
                        plot(ts,PSTH(chan,:),'k'); xlabel('Time (s)'); ylabel('Spike rate (Hz)')
                        try
                            set(sp,'YDir','normal')
                            set(sp,'xtick',freq_tick,...
                                'xticklabel',freq_val,...
                                'FontSize',7,...
                                'XTickLabelRotation',45)
                            set(pst,'FontSize',7)
                        catch
                        end
                        
                        chan  = chan + 1;
                        if chan > nChan
                            break
                        end
                    end
                    if chan > nChan
                        break
                    end
                end
                if chan > nChan
                    break
                end
            end
            % Plot an overall depth fig
            for s = unique(sitedepths(:,3))'
                figure;
                % Characteristic freq
                ax1 = subplot(1,4,1);                
                for i = find(sitedepths(:,3)==s)'
                    plot(storeQ30(i,:),repmat(sitedepths(i,2),1,2),'-','color',[33, 131, 128]/255,'LineWidth',2); hold on
                    plot(storeQ10(i,:),repmat(sitedepths(i,2),1,2),'-','color',[216, 17, 89]/255,'LineWidth',3); hold on
                    if ~isnan(OnsetLatencys(i))
                         plot(storeCF(i),sitedepths(i,2),'kd','markerfacecolor','k','MarkerSize',5); xlabel('Freq'); hold on
                    else
                        plot(storeCF(i),sitedepths(i,2),'.','color',[155, 161, 163]/255,'MarkerSize',15); xlabel('Freq'); hold on
                    end
                end
                legend({'Q30','Q10'})
                set(gca,'XScale','log','xtick',round(newfreqs(freq_tick),2).*1e3,'ytick',min(sitedepths(sitedepths(:,3)==s,2)):250:max(sitedepths(sitedepths(:,3)==s,2)))
                ylabel('Depth from tip (um)')

                % Onset latencys and lengths
                ax2 = subplot(1,4,2);
                plot(OnsetLatencys(find(sitedepths(:,3)==s)),sitedepths(find(sitedepths(:,3)==s),2),'kd','markerfacecolor','k','MarkerSize',5); xlabel('Onset latency (s)');yticks(min(sitedepths(:,2)):250:max(sitedepths(:,2)))
                ax3 = subplot(1,4,3);
                plot(OnsetLengths(find(sitedepths(:,3)==s)),sitedepths(find(sitedepths(:,3)==s),2),'kd','markerfacecolor','k','MarkerSize',5); xlabel('Onset lengths (s)');yticks(min(sitedepths(:,2)):250:max(sitedepths(:,2)))
                linkaxes([ax1 ax2 ax3],'y')
                ax4 = subplot(1,4,4);
                imagesc(flipud(storeBounds(sitedepths(:,3)==s & ~isnan(OnsetLatencys)',:)*10)); a = colorbar; a.Label.String = 'Level (dB)';
                colormap(flipud(bone))
                set(gca,'ytick',1:size(storeBounds(sitedepths(:,3)==s & ~isnan(OnsetLatencys)',:),1),...
                    'yticklabel',flipud(sitedepths(sitedepths(:,3)==s & ~isnan(OnsetLatencys)',2)),...
                    'xtick',freq_tick,...
                    'xticklabel',freq_val)
                xlabel('Freq')
                sgtitle(['Shank ',num2str(s)])
            end
        end

        
    end
end

