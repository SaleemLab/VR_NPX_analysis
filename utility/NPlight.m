classdef NPlight
    %NPLIGHT Code to view data from the light levels

    methods (Static)
       function CSD(binName,binPath)
            ISI = 0.4;
            window = 4.5;
            % Get all the info
            [trialTimes,~] = NPadmin.getTrialTimes(binName,binPath);
            [LFP,config,info] = getNPLFP(binName,binPath,[]);
            LFPStartSamples = floor((trialTimes-ISI)*info.fs);
            LFPEndSamples = floor((ISI+window)*info.fs);
            stimTable = StimRef.getTrialInfo(binName, binPath,trialTimes);
            
            for currPL = unique(stimTable(:,4))'
                for side = 0:1
%                     subplot(1,2,side+1)
                    currIdxs = find(stimTable(:,4)==currPL & stimTable(:,8)==side);
                    % Create a matrix where dim1 = samples. dim2 = channels, dim3 = trials
                    LFPmat = [];
                    for currTrial = currIdxs'
                        LFPmat(:,:,currTrial == currIdxs) = LFP(LFPStartSamples(currTrial):LFPStartSamples(currTrial)+LFPEndSamples,:);
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
                    sgtitle(['PL: ',num2str(currPL),' REG: ',num2str(side)])
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
            end
        end
    end
end

