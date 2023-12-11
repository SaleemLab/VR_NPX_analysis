classdef NPksAdmin
    %NPKSADMIN admin functions to extract spikes from kilosort/phy
    %   Detailed explanation goes here
    
    methods (Static)
        % Script to run kilosort on everything that hasn't been sorted yet
        function batchKilosort
            
            answer = questdlg('Do you want to overwrite previous spike sorting?',...
                'Yes','No');
            switch answer
                case 'Yes'
                    overwrite = 1;
                case 'No'
                    overwrite = 0;
                case 'Cancel'
                    disp('Cancelled spikesorting')
                    return
            end
            
            % Add paths
            pathToYourConfigFile = '.\toolboxes\Kilosort-main\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
            inDataDir = 'E:\NeuralData\Neuropixel Data\';
            folders = dir([inDataDir,'F*']);
            for currFolder = 1:size(folders,1)
                currFerret = folders(currFolder).name;
                
                currInDataDir = [inDataDir,currFerret,'\'];
                
                NpFolders = dir([currInDataDir,'*_g*']);
                
                for currNPFolder = 1:size(NpFolders,1)
                    npDir = [NpFolders(currNPFolder).name,'\'];
                    IMECdirs = dir([currInDataDir,npDir,[npDir(1:end-1),'_imec*']]);
                    for currIMEC = 1:size(IMECdirs,1)
                        currDir = [currInDataDir,npDir,IMECdirs(currIMEC).name];
                        
                        % Do some checks so we arn't respikesorting
                        %Firstly check if the file is 0 (means it crashed)
                        tempFilesSize = dir([currDir,'\*.ap.bin']);
                        
                        % Delete kilosort folder if overwriting
                        if overwrite == 1 && isfolder([currDir,'\kilosort3'])
                                dos_cmd = sprintf( 'rmdir /S /Q "%s"', [currDir,'\kilosort3']);
                                system(dos_cmd);
                        end
                        tempFileCheck = dir([currDir,'\kilosort3\params.py']);
                        recFileCheck = dir([currInDataDir,npDir(1:end-1),'\*_t0.nidq.bin']);
                                                try
                        if tempFilesSize.bytes ~= 0 && (isempty(tempFileCheck) || overwrite == 1) && ~isempty(recFileCheck)
                            % Get directory
                            rootZ = currDir; % the raw data binary file is in this folder
                            %                               % path to temporary binary file (same size as data, should be on fast SSD)
                            rootH = rootZ; % Used to be on d drive, but not temp as phy needs it to read in the waveforms
                            
                            % check if theres a temp file and if so
                            % delete it (to avoid wierd waveforms in
                            % phy)
                            try
                                delete(fullfile(rootH,'temp_wh.dat'))
                            end
                            
                            binName = [npDir(1:end-1),'_t0.nidq.bin'];
                            binPath = [currInDataDir,npDir(1:end-1)];
                            filename = [currInDataDir,npDir,binName];
                            meta     = NPadmin.ReadMeta(binName, binPath);
                            APbinDir = dir([currDir,'/*_t0.imec*.ap.bin']);
                            APbinName = APbinDir.name;
                            metaAP     = NPadmin.ReadMeta(APbinName, currDir);
                            
                            % Need to find a way to check the channel map used
                            [~,imroFileName] = fileparts(metaAP.imRoFile);
                            chanMapFile = [imroFileName,'.mat'];
                            ops.trange = [0 Inf]; % time range to sort
                            ops.NchanTOT    = 385; % total number of channels in your recording % Was 384
                            
                            
                            run(fullfile(pathToYourConfigFile, 'configFile384.m'))
                            ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
                            ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);
                            
                            
                            
                            fprintf('Looking for data inside %s \n', rootZ)
                            
                            % main parameter changes from Kilosort2 to v2.5
                            ops.sig        = 20;  % spatial smoothness constant for registration
                            ops.fshigh     = 300; % high-pass more aggresively
                            ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
                            
                            % main parameter changes from Kilosort2.5 to v3.0
                            ops.Th       = [9 9];
                            
                            % is there a channel map file in this folder?
                            fs = dir(fullfile(rootZ, 'chan*.mat'));
                            if ~isempty(fs)
                                ops.chanMap = fullfile(rootZ, fs(1).name);
                            end
                            % find the binary file
                            fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
                            ops.fbinary = fullfile(rootZ, fs(1).name);
                            
                            rez                = preprocessDataSub(ops);
                            rez                = datashift2(rez, 1);
                            
                            [rez, st3, tF]     = extract_spikes(rez);
                            
                            rez                = template_learning(rez, tF, st3);
                            
                            [rez, st3, tF]     = trackAndSort(rez);
                            
                            rez                = final_clustering(rez, tF, st3);
                            
                            rez                = find_merges(rez, 1);
                            
                            rootZ = fullfile(rootZ, 'kilosort3');
                            mkdir(rootZ)
                            rezToPhy2(rez, rootZ);
                            disp('Successfully kilosorted')
                            
                            % Clear a bunch of variables
                            clear tf st3 rez ops meta metaAP
                        end
                                                catch
%                                                     kaja = [];
                                                    fprintf('Failed at %s \n', rootZ)
                        %
                                                end
                    end
                end
            end
        end
        % Getting spikes
        function [spike] = extractPhySpikes(currDir)
            % Get out spike data
            sp = NPksAdmin.loadKSdir(currDir);
            disp(currDir)
            spikeTimes = sp.st;
            spikeTemps = sp.spikeTemplates;

            
            
            % Due to kilsorot 3 not having a pc features output file need
            % to fudge around alot
            [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
                NPksAdmin.templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
            
            [~,max_site] = max(max(abs(tempsUnW),[],2),[],3); % the maximal site for each template
            spikeSites = max_site(spikeTemps+1); % Index from 0
            
            % Add in cluster groups in the second column
            % - 0 = noise
            % - 1 = mua
            % - 2 = good
            % - 3 = unsorted
            % cids is length nClusters, the cluster ID numbers
            % cgs is length nClusters, the "cluster group":
            for i = 1:size(spikeTemps,1)
                try
                    spikeTemps(i,2) = sp.cgs(sp.cids == spikeTemps(i,1));
                catch % If it can't find it then its been deleted and was classed as noise
                    spikeTemps(i,2) = 0;
                end
            end
            
            spike.Times = spikeTimes;
            spike.Amps = spikeAmps;
            spike.Depths = spikeDepths; % This is the depth computed as a centre of mass (no drift)
            spike.Sites = spikeSites;
            spike.Ycoords = sp.ycoords;
            spike.Xcoords = sp.xcoords;
            spike.templateYpos = templateYpos;
            spike.tempAmps = tempAmps;
            spike.tempPeakWF = tempPeakWF;
            spike.templates = spikeTemps;
        end
        
        % Phy helpers taken from Cortex lab
        % (https://github.com/cortex-lab/spikes/)
        function spikeStruct = loadKSdir(ksDir, varargin)
            
            if ~isempty(varargin)
                params = varargin{1};
            else
                params = [];
            end
            
            if ~isfield(params, 'excludeNoise')
                params.excludeNoise = true;
            end
            if ~isfield(params, 'loadPCs')
                params.loadPCs = false;
            end
            
            % load spike data
            
            spikeStruct = NPksAdmin.loadParamsPy(fullfile(ksDir, 'params.py'));
            
            ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
            st = double(ss)/spikeStruct.sample_rate;
            spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
            
            if exist(fullfile(ksDir, 'spike_clusters.npy'))
                clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
            else
                clu = spikeTemplates;
            end
            
            tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));
            
            if params.loadPCs
                pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
                pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
            else
                pcFeat = [];
                pcFeatInd = [];
            end
            
            cgsFile = '';
            if exist(fullfile(ksDir, 'cluster_groups.csv'))
                cgsFile = fullfile(ksDir, 'cluster_groups.csv');
            end
            if exist(fullfile(ksDir, 'cluster_group.tsv'))
                cgsFile = fullfile(ksDir, 'cluster_group.tsv');
            end
            if ~isempty(cgsFile)
                [cids, cgs] = NPksAdmin.readClusterGroupsCSV(cgsFile);
                
                if params.excludeNoise
                    noiseClusters = cids(cgs==0);
                    
                    st = st(~ismember(clu, noiseClusters));
                    spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
                    tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
                    
                    if params.loadPCs
                        pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
                        %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
                    end
                    
                    clu = clu(~ismember(clu, noiseClusters));
                    cgs = cgs(~ismember(cids, noiseClusters));
                    cids = cids(~ismember(cids, noiseClusters));
                    
                    
                end
                
            else
                clu = spikeTemplates;
                
                cids = unique(spikeTemplates);
                cgs = 3*ones(size(cids));
            end
            
            
            coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
            ycoords = coords(:,2); xcoords = coords(:,1);
            temps = readNPY(fullfile(ksDir, 'templates.npy'));
            
            winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));
            
            spikeStruct.st = st;
            spikeStruct.spikeTemplates = spikeTemplates;
            spikeStruct.clu = clu;
            spikeStruct.tempScalingAmps = tempScalingAmps;
            spikeStruct.cgs = cgs;
            spikeStruct.cids = cids;
            spikeStruct.xcoords = xcoords;
            spikeStruct.ycoords = ycoords;
            spikeStruct.temps = temps;
            spikeStruct.winv = winv;
            spikeStruct.pcFeat = pcFeat;
            spikeStruct.pcFeatInd = pcFeatInd;
        end
        
        function [cids, cgs] = readClusterGroupsCSV(filename)
            %function [cids, cgs] = readClusterGroupsCSV(filename)
            % cids is length nClusters, the cluster ID numbers
            % cgs is length nClusters, the "cluster group":
            % - 0 = noise
            % - 1 = mua
            % - 2 = good
            % - 3 = unsorted
            
            fid = fopen(filename);
            C = textscan(fid, '%s%s');
            fclose(fid);
            
            cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
            ise = cellfun(@isempty, cids);
            cids = [cids{~ise}];
            
            isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
            isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
            isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
            cgs = zeros(size(cids));
            
            cgs(isMUA) = 1;
            cgs(isGood) = 2;
            cgs(isUns) = 3;
        end
        
        % Mostly for just getting the sp structure
        function [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksDir, localizedSpikesOnly)
            
            if nargin < 2
                localizedSpikesOnly = false;
            end
            
            clear params
            params.loadPCs = true;
            sp = loadKSdir(ksDir, params);
            
            if localizedSpikesOnly % go over all templates and check which are not localized (in space)
                localizedTemplates = false(size(sp.temps,1), 1);
                for t = 1:size(sp.temps,1)
                    M = max(max(abs(squeeze(sp.temps(t,:,:)))));
                    ch = find(max(abs(squeeze(sp.temps(t,:,:)))) > 0.5*M); % the channels where the template has significant weight
                    localizedTemplates(t) = max(ch) - min(ch) <= 20; % all channels at most 20 apart?
                    
                    %assert(sum(sum(squeeze(sp.temps(t,:,:)).^2)) == 1, 'template norm not 1?!')
                    %com = sp.ycoords .* max(abs(squeeze(sp.temps(t,:,:)))).^2 / sum(max(abs(squeeze(sp.temps(t,:,:)))).^2)
                end
                localizedTemplates = uint32(find(localizedTemplates) - 1); % the numbers of localized templates (indexing starts from 0)
                i = ismember(sp.spikeTemplates, localizedTemplates);
                sp.st = sp.st(i);
                sp.spikeTemplates = sp.spikeTemplates(i);
                sp.clu = sp.clu(i);
                sp.tempScalingAmps = sp.tempScalingAmps(i);
                sp.pcFeat = sp.pcFeat(i,:,:);
                clear i localizedTemplates M ch t
            end
            
            ycoords = sp.ycoords;
            pcFeat = sp.pcFeat;
            pcFeat = squeeze(pcFeat(:,1,:)); % take first PC only
            pcFeat(pcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there.
            pcFeatInd = sp.pcFeatInd;
            spikeTemps = sp.spikeTemplates;
            
            temps = sp.temps;
            winv = sp.winv;
            tempScalingAmps = sp.tempScalingAmps;
            spikeTimes = sp.st;
            
            % compute center of mass of these features
            
            % which channels for each spike?
            spikeFeatInd = pcFeatInd(spikeTemps+1,:);
            
            % ycoords of those channels?
            spikeFeatYcoords = ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12
            
            % center of mass is sum(coords.*features)/sum(features)
            spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);
            
            
            % for plotting, we need the amplitude of each spike, both so we can draw a
            % threshold and exclude the low-amp noisy ones, and so we can color the
            % points by the amplitude
            
            [spikeAmps, ~, templateYpos, tempAmps, tempsUnW] = ...
                templatePositionsAmplitudes(temps, winv, ycoords, spikeTemps, tempScalingAmps);
            
            %[~,max_site] = max(max(abs(temps),[],2),[],3); % the maximal site for each template
            % one could potentially use the unwhitened templates, but that shouldn't really change the results
            [~,max_site] = max(max(abs(tempsUnW),[],2),[],3); % the maximal site for each template
            spikeSites = max_site(spikeTemps+1);
            
            if isfield(sp, 'gain') && ~isempty(sp.gain) % could put this field in your params.py
                % spikeAmps = spikeAmps*0.6/512/500*1e6; % convert to uV
                spikeAmps = spikeAmps*sp.gain; % convert to uV
            end
            
            if localizedSpikesOnly  % above we already removed non-localized templates, but that on its own is insufficient
                b = regress(spikeSites, spikeDepths); % for IMEC probe adding a constant term kills the regression making the regressors rank deficient
                i = abs(spikeSites - b*spikeDepths) <= 5;
                spikeTimes  = spikeTimes(i);
                spikeAmps   = spikeAmps(i);
                spikeDepths = spikeDepths(i);
                spikeSites  = spikeSites(i);
            end
            spikeSites = uint16(spikeSites);
        end
        
        function S = loadParamsPy(fn)
            % Loads a phy-style "params.py" into a matlab struct. The params.py is
            % python code but just a series of assignments, so most of it will run
            % directly in matlab
            %
            % now based on text2struct from James Jun, 2016-12-31, modified slightly by
            % N. Steinmetz.
            
            fid = fopen(fn, 'r');
            mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
            fclose(fid);
            csName = mcFileMeta{1};
            csValue = mcFileMeta{2};
            S = struct();
            for i=1:numel(csName)
                vcName1 = csName{i};
                if vcName1(1) == '~', vcName1(1) = []; end
                try
                    if csValue{i}(1)==''''; % if first character of the value is a single quote
                        % then it's a string and we can evaluate without adding a
                        % single quote
                        eval(sprintf('%s = %s;', vcName1, csValue{i}));
                    else
                        eval(sprintf('%s = ''%s'';', vcName1, csValue{i}));
                    end
                    eval(sprintf('num = str2double(%s);', vcName1));
                    if ~isnan(num)
                        eval(sprintf('%s = num;', vcName1));
                    end
                    eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
                catch
                    fprintf('%s = %s error\n', csName{i}, csValue{i});
                end
            end
            
            % handle special case of true/false
            fnames = fieldnames(S);
            for f = 1:length(fnames)
                if strcmp(S.(fnames{f}), 'True') || strcmp(S.(fnames{f}), 'true')
                    S.(fnames{f})=true;
                elseif strcmp(S.(fnames{f}), 'False') || strcmp(S.(fnames{f}), 'false')
                    S.(fnames{f})=false;
                end
            end
            
            
        end
        
        function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
            % function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
            %
            % Compute some basic things about spikes and templates
            %
            % outputs:
            % - spikeAmps is length nSpikes vector with amplitude in unwhitened space
            % of every spike
            % - spikeDepths is the position along the probe of every spike (according
            % to the position of the template it was extracted with)
            % - templateDepths is the position along the probe of every template
            % - templateAmps is the amplitude of each template
            % - tempsUnW are the unwhitened templates
            % - templateDuration is the trough-to-peak time (in samples)
            % - waveforms: returns the waveform from the max-amplitude channel
            %
            % inputs:
            % - temps, the templates (nTemplates x nTimePoints x nChannels)
            % - winv, the whitening matrix (nCh x nCh)
            % - ycoords, the coordinates of the channels (nCh x 1)
            % - spikeTemplates, which template each spike came from (nSpikes x 1)
            % - tempScalingAmps, the amount by which the template was scaled to extract
            % each spike (nSpikes x 1)
            
            % unwhiten all the templates
            tempsUnW = zeros(size(temps));
            for t = 1:size(temps,1)
                tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
            end
            
            % compute the biggest absolute value within each template (obsolete)
            % absTemps = abs(tempsUnW);
            % tempAmps = max(max(absTemps,[],3),[],2);
            
            % The amplitude on each channel is the positive peak minus the negative
            tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));
            
            % The template amplitude is the amplitude of its largest channel (but see
            % below for true tempAmps)
            tempAmpsUnscaled = max(tempChanAmps,[],2);
            
            % need to zero-out the potentially-many low values on distant channels ...
            threshVals = tempAmpsUnscaled*0.3;
            tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;
            
            % ... in order to compute the depth as a center of mass
            templateDepths = sum(bsxfun(@times,tempChanAmps,ycoords'),2)./sum(tempChanAmps,2);
            
            % assign all spikes the amplitude of their template multiplied by their
            % scaling amplitudes (templates are zero-indexed)
            spikeAmps = tempAmpsUnscaled(spikeTemplates+1).*tempScalingAmps;
            
            % take the average of all spike amps to get actual template amps (since
            % tempScalingAmps are equal mean for all templates)
            ta = clusterAverage(spikeTemplates+1, spikeAmps);
            tids = unique(spikeTemplates);
            tempAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike
            tempAmps = tempAmps'; % for consistency, make first dimension template number
            
            % Each spike's depth is the depth of its template
            spikeDepths = templateDepths(spikeTemplates+1);
            
            % Get channel with largest amplitude, take that as the waveform
            [~,max_site] = max(max(abs(temps),[],2),[],3);
            templates_max = nan(size(temps,1),size(temps,2));
            for curr_template = 1:size(temps,1)
                templates_max(curr_template,:) = ...
                    temps(curr_template,:,max_site(curr_template));
            end
            waveforms = templates_max;
            
            % Get trough-to-peak time for each template
            [~,waveform_trough] = min(templates_max,[],2);
            [~,templateDuration] = arrayfun(@(x) ...
                max(templates_max(x,waveform_trough(x):end),[],2), ...
                transpose(1:size(templates_max,1)));
            
            
            
            
            
        end
    end
end
