classdef NPoverview
    %NPOVERVIEW Overview functions
    %
    properties
    end
    
    methods (Static)
        
        function LFPoverview(binName,binPath)
            [LFP,config,info] = getNPLFP(binName,binPath,240);
            
            nfft = 2^(nextpow2(info.fs));
            win  = hanning(nfft);
            [pxx,fxx] = pwelch(LFP,win,[],nfft,info.fs);
            
            freqs  = [3 7;30 70;150 250;300 350;550 700];
            
            [P, F] = meanPSDfreqrange(pxx,fxx,freqs);
            normP = P ./ max(P);
            
            
            % Do some plotting
            totalShanks = numel(unique(config.Shank));
            selectedShanks = unique(config.Shank)';
            figure;
            tStart = 10;
            tEnd = tStart+5;
            for s = 1:totalShanks
                ax1 = subplot(1,totalShanks*3,s);
                
                currChans = config.Channel(config.Shank==selectedShanks(s));
                allLFProwIdx = find(ismember(config.Channel,currChans));
                allElect_Sel = config.Electrode(allLFProwIdx);
                currChans = currChans(1:5-totalShanks:numel(currChans)); % Choose a subset for plotting
                LFProwIdx = find(ismember(config.Channel,currChans));
                elec_selection = config.Electrode(LFProwIdx);
                
                t = tStart : 1/info.fs :tEnd;  % rough timeline
                NPplot.plotTraceAcrossChans(gca,LFP(floor(info.fs*tStart):floor(info.fs*tStart)+(numel(t)-1),LFProwIdx),'xvec',t,'yvec',elec_selection/2,'scalefactor',0.01)
                ylabel('electrode')
                xlabel('time (s)')
                if info.probeType == 24
%                     set(ax1, 'ylim',[-10,650], 'ytick',0:50:650,'yticklabel',0:100:1300,'xlim',[t(1) t(end)])
%                     electrodes = 1:1280; % all electrodes
                              set(ax1, 'ylim',[-10,650], 'ytick',0:50:650,'yticklabel',0:100:1300,'xlim',[t(1) t(end)])
                    electrodes = 1:1280; % all electrodes
                else
%                     set(ax1, 'ylim',[-10,500], 'ytick',0:50:450,'yticklabel',0:100:900,'xlim',[t(1) t(end)])
%                     electrodes = 1:960; % all electrodes
                    set(ax1, 'ylim',[-10,200], 'ytick',0:50:450,'yticklabel',0:100:400,'xlim',[t(1) t(end)])
                    electrodes = 1:960; % all electrodes
                    
                end
                title(['LFP Shank ',num2str(selectedShanks(s))])
                
                ax2 = subplot(1,totalShanks*3,s+totalShanks);
                
                colsG = grayscale_based_on_num_input(size(freqs,1));
                colsG(:,1) = flip(colsG(:,1));
                
                new_nP = NaN(max(electrodes),size(normP,2));
                new_nP(ismember(electrodes,elec_selection),:)=normP(LFProwIdx,:);
                for n = 1:size(freqs,1)
                    txt = [num2str(round(F(n,1))) '-' num2str(round(F(n,2))) 'Hz'];
                    plot(new_nP(ismember(electrodes,elec_selection),n),elec_selection/2,'Marker','.','color',colsG(n,:),'LineWidth',2,'DisplayName',txt);hold on
                end
                if info.probeType == 24
                    set(ax2,'ylim',[-10 650],'ycolor','none','xtick',[0 1])
                else
%                     set(ax2,'ylim',[-10 500],'ycolor','none','xtick',[0 1])
                    set(ax2,'ylim',[-10 200],'ycolor','none','xtick',[0 1])
                end
                xlabel('Norm Power')
                legend show
                title(['Norm P Shank ',num2str(selectedShanks(s))])
                
                % Plot theta phase shift
                % get phase shift
                refChan = 1;
                %% Masa: maybe need to change phase shift?
                [phaseDiff,phaseDiffStd] = getPhaseShift(LFP(:,allLFProwIdx),3,7,info.fs,refChan);
                ax3 = subplot(1,totalShanks*3,s+(totalShanks*2));
                plot(phaseDiff,allElect_Sel/2,'Marker','.','color',colsG(1,:),'LineWidth',2);hold on
                patch([phaseDiff-phaseDiffStd flip(phaseDiff+phaseDiffStd)],[allElect_Sel/2;flip(allElect_Sel/2)],colsG(1,:),'EdgeAlpha',0,'FaceAlpha',.5)
                xlabel(['Phase shift (degrees) from elec ',num2str(allElect_Sel(refChan))])
                if info.probeType == 24
                    set(ax3,'ylim',[-10 650],'ycolor','none')
                else
                    set(ax3,'ylim',[-10 200],'ycolor','none')
%                     set(ax3,'ylim',[-10 500],'ycolor','none')
                end
                title(['Phase shift Shank ',num2str(selectedShanks(s))])
                
            end
        end
        function LFPcorrmat(binName,binPath)
            [LFP,config,info] = getNPLFP(binName,binPath,120);
            r = corrcoef(LFP,'Rows','complete');
            sites = config.Electrode;
            
            figure; h = pcolor(r);
            set(gca,'clim',[-1 1],'fontsize', 8)
            set(h,'EdgeColor','none')
            axis square;
            xticks(1:10:size(r,1));yticks(1:10:size(r,1))
            xticklabels(num2str(sites(1:10:size(r,1))));xtickangle(270)
            yticklabels(num2str(sites(1:10:size(r,1))));
            title('Correlation matrix'); xlabel('Electrode');ylabel('Electrode')
            % Add shank markers
            if info.probeType == 24
                shanks = config.Shank;
                shankIdx = find(diff(shanks));
                for i = 1:numel(shankIdx)
                    xline(shankIdx(i),'k')
                    yline(shankIdx(i),'k')
                end
                figure;
                % Plot each shank individually as well
                totalShanks = 4;
                for i = 1:totalShanks
                    subplot(2,2,i)
                    currChans = find(config.Shank==i);
                    currR = r(currChans,currChans);
                    h = pcolor(currR);
                    sites = config.Electrode(currChans);
                    set(h,'EdgeColor','none')
                    set(gca,'clim',[-1 1],'fontsize', 8)
                    axis square;
                    xticks(1:10:size(currR,1));yticks(1:10:size(currR,1))
                    xticklabels(num2str(sites(1:10:size(currR,1))));xtickangle(270)
                    yticklabels(num2str(sites(1:10:size(currR,1))));
                    xlabel('Electrode');ylabel('Electrode')
                    title(['Shank ',num2str(i)])
                end
            end
        end
        
        function LFPonset(binName,binPath)
            window = 0.6;
            ISI = 0.2;
            
            [trialTimes,~] = NPadmin.getTrialTimes(binName,binPath);
            
            [LFP,config,info] = getNPLFP(binName,binPath,[]);
            LFPStartSamples = floor((trialTimes-ISI)*info.fs);
            LFPEndSamples = floor((ISI+window)*info.fs);
            
            % Create a matrix where dim1 = samples. dim2 = channels, dim3 = trials
            LFPmat = [];
            for currTrial = 1:numel(LFPStartSamples)
                LFPmat(:,:,currTrial) = LFP(LFPStartSamples(currTrial):LFPStartSamples(currTrial)+LFPEndSamples,:);
            end
            
            LFPmean = mean(LFPmat,3);
            scaleFactor = 0.5;
            figure;
            t = -ISI:1/info.fs:window;
            if info.probeType == 24
                for s = 1:4
                    LFPidx = find(config.Shank==s);
                    subplot(1,4,s);
                    NPplot.plotTraceAcrossChans(gca,LFPmean(:,LFPidx),'xvec',t, 'yvec',config.Electrode(LFPidx), 'scalefactor',scaleFactor)
                    ylim([-10 1300]); ylabel('Electrode'); xlabel('Time (s)');
                    title(['Shank ',num2str(s)])
                end
            else
                NPplot.plotTraceAcrossChans(gca,LFPmean(:,1:2:end),'xvec',t,'yvec', config.Electrode(1:2:end), 'scalefactor',scaleFactor)
                ylim([-10 400]); ylabel('Electrode'); xlabel('Time (s)'); xlim([-ISI window])
%                 ylim([-10 950]); ylabel('Electrode'); xlabel('Time (s)'); xlim([-ISI window])
                title('Average LFP on trial Onset')
            end
            
        end
        
        function LFPonsetCSD(binName,binPath)
            window = 0.4;
            ISI = 0.2;
            
            [trialTimes,~] = NPadmin.getTrialTimes(binName,binPath);
            
            [LFP,config,info] = getNPLFP(binName,binPath,[]);
            LFPStartSamples = floor((trialTimes-ISI)*info.fs);
            LFPEndSamples = floor((ISI+window)*info.fs);
            
            % Create a matrix where dim1 = samples. dim2 = channels, dim3 = trials
            LFPmat = [];
            for currTrial = 1:numel(LFPStartSamples)
                LFPmat(:,:,currTrial) = LFP(LFPStartSamples(currTrial):LFPStartSamples(currTrial)+LFPEndSamples,:);
            end
            
            LFPmean = mean(LFPmat,3);
            
            % Baseline correct
            for i = 1:size(LFPmean,2)
                baseline = mean(LFPmean(1:floor(ISI*info.fs),i));
                LFPmean(:,i) = LFPmean(:,i)-baseline;
            end
            
            t = -ISI:1/info.fs:window;
            
            totalShanks = numel(unique(config.Shank));
            selectedShanks = unique(config.Shank)';
            figure
            for s = 1:totalShanks
                subplot(1,totalShanks,s)
                xval = mode(config.Ks_xcoord(config.Shank == selectedShanks(s)));
                LFPselect = find(config.Ks_xcoord==xval & config.Shank == selectedShanks(s));
                ycoords = config.Ks_ycoord(LFPselect);
                spacing = round(mean(diff(ycoords)));
                [CSDoutput]=LFPanalysis.CSD(LFPmean(:,LFPselect),info.fs,spacing);
                contourf(t,ycoords,CSDoutput',40,'linecolor','none'); hold on
                xlabel('Time (s)'); ylabel('Electrode')
                colorbar
                cBmap = cbrewer('div','RdBu',5);
                cBmapInt = interpolate_cbrewer(cBmap, 'linear', 20);
                colormap(redblue);
                Clim = max([max(CSDoutput) abs(min(CSDoutput))]);
                caxis([-Clim Clim]);
                title(['Shank: ',num2str(selectedShanks(s)),'. Source = blue, Sink = red'])
            end  
        end
        
        function Spikeonset(binName,binPath)
            % check if its been sorted by Phy
            sortedCheck = dir([binPath,'kilosort3\','phy.log']);
            if isempty(sortedCheck)
                errordlg('Not sorted yet')
                return
            end
            window = 0.6;
            ISI = 0.2;
            
            [StartTimes,~] = NPadmin.getTrialTimes(binName,binPath);
            spike = NPksAdmin.extractPhySpikes([binPath,'kilosort3\']);
            meta = NPadmin.ReadMeta(binName, binPath);
            EndTimes = ISI+window;
            
            % Associate spike times with spike Sites
            spikeData(:,1) = spike.Times;
            spikeData(:,2) = spike.Sites;
            % Associate with depth
            spikeData(:,3) = spike.Ycoords(spike.Sites);
            % Associate with shank
            if str2double(meta.imDatPrb_type) == 24
                % Need to figure out shank from xcoords
                shanks = NPadmin.xcoordsToShank(spike.Xcoords);
                spikeData(:,4) = shanks(spike.Sites);
            else
                spikeData(:,4) = ones(size(spikeData,1),1);
            end
            
            disp('Getting spike data...')
            for currShank = 1:max(unique(spikeData(:,4)))
                spikeDepths = (unique(spikeData(spikeData(:,4)==currShank,3))/10)'; % Dividing by ten to make it easier for the array
                if ~isempty(spikeDepths)
                    outArray = [];ytickNames = [];selectLocations = []; outTaso = [];tasoCount = [];
                    figure('name',['Shank: ',num2str(currShank)]);
                    for currDepth = spikeDepths
                        Sites = unique(spikeData(spikeData(:,3)==(currDepth*10)&spikeData(:,4)==currShank,2));
                        for currSite = Sites'
                            currSpikeTimes = spikeData(spikeData(:,3)==(currDepth*10)&spikeData(:,2)==currSite&spikeData(:,4)==currShank,1);
                            
                            [y,ts,taso]= NPplot.NPgetPSTH(currSpikeTimes,StartTimes,ISI,window,0.005);
                            % Baseline on the ISI
                            y = y-nanmean(y(ts<ISI));
                            %         outArray(currDepth==spikeDepths,:) = y;
                            y = zscore(y);
                            outArray = [outArray;y];
                            outTaso = [outTaso;taso];
                            % count to know how many trials for each depth
                            tasoCount = [tasoCount size(outTaso,1)];
                            ytickNames = [ytickNames currDepth];
                            selectLocations = [selectLocations; spike.Xcoords(currSite) spike.Ycoords(currSite)];
                        end
                    end
                    
                    xtickBinWidth = 0.1;
                    ax1 = subplot(1,3,1);
                    disp('Plotting...')
                    imagesc(flipud(outArray)); title('Depth = ycoord')
                    displayTimes = round([-ISI:xtickBinWidth:window],2);
                    displayDepths = spikeDepths(1:0.1*numel(spikeDepths):end);
                    xtickIdx = find(ismember(round(ts*10,2),displayTimes*10));
                    xticks(xtickIdx); xticklabels(ts(xtickIdx));
                    % yticks(ytickIdx);yticklabels(spikeDepths(ytickIdx)*10)
                    
                    yticks([1:5:size(outArray,1)]);yticklabels(flip(ytickNames(1:5:end)*10))
                    xlabel('Time (s)'); ylabel('Depth from tip (um)'); h = colorbar; ylabel(h, 'Z-scored Spike Rate (Hz)')
                    ax2 = subplot(1,3,2);
                    % Plot the line PSTHs
                    for i = 1:size(outArray,1)
                        % Scale factor
                        scaleFac = 10;
                        plot(ts,(outArray(i,:)/scaleFac)+i,'k-'); hold on
                    end
                    %
                    title('PSTH')
                    yticks([1:5:size(outArray,1)]);yticklabels(ytickNames(1:5:end)*10)
                    xlabel('Time (s)'); ylabel('Depth from tip (um)'); ylabel(h, 'Z-scored Spike Rate (Hz)')
                    ax3 = subplot(1,3,3);
                    outTaso(outTaso==0)=NaN;
                    disp('Plotting raster, may take some time')
                    % Plot the rasters
                    for i = 1:size(outTaso,1)
                        plot(outTaso(i,:)+i,'k.','MarkerSize',1); hold on
                    end
                    yticks([tasoCount(1:5:end)]);
                    yticklabels(ytickNames(1:5:end)*10); ylim([0 max(tasoCount)])
                    title('Raster'); xlabel('Time (s)'); ylabel('Depth from tip (um)')
                    xtickIdx = find(ismember(round(([0:size(outTaso,2)]/30000)-ISI,5),displayTimes));
                    xticks(xtickIdx);
                    xticklabels(displayTimes);
                    sgtitle(['Shank ',num2str(currShank)])
                end
            end
        end
        
    end
end

