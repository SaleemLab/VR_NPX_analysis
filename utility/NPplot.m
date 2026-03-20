classdef NPplot
    %NPPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        % Getting functions
        function [y,ts,taso,spikeRate]= NPgetPSTH(spikeTimes,StartSamples,ISI,stimLength,binWidth)
            fs = 30000;
            StartSamples = floor(StartSamples*fs); % To convert the time into samples
            spikeTimes = floor(spikeTimes*fs);
            
            BWfs = floor(binWidth*fs);
            nBins = floor(((stimLength + ISI)*fs)/BWfs);
            taso = zeros(1,ceil(((stimLength + ISI)*fs)));
            
            
            for i = 1:length(StartSamples)
                begSample = StartSamples(i) - ceil(ISI*fs);
                endSample = begSample + ceil((stimLength + ISI)*fs);
                currSpikeTimes = spikeTimes(spikeTimes > begSample & spikeTimes < endSample,1);
                % convert to time
                currSpikeTimes = currSpikeTimes - begSample;
                taso(i,currSpikeTimes) = ones(1,length(currSpikeTimes));
                
                % Then figure out some firing rates
                for currBin = 1:nBins
                    spikeRate(i,currBin) = sum(taso(i,(currBin*BWfs)-BWfs+1:currBin*BWfs))./binWidth;
                end
                
            end
            
            y = mean(spikeRate);
            ts = -ISI:binWidth:((length(y)-1)*binWidth)-ISI; % Old version
            
        end
        % Plotting functions
        function scalefactor = plotTraceAcrossChans(ax1, signal, varargin)
            %
            % scalefactor = plot_trace_across_chans(ax1, signal, varargin)
            % name pair inputs: xvec,yvec,zval,,scalefactor,color,linewidth, linestyle
            % Soraya Dunn 2020
            %
            
            if ndims(signal) > 2 %#ok<*ISMAT> % supress suggestion
                signal = squeeze(signal);
                if ndims(signal) > 2
                    disp('cannot use plot_trace_across_chans.m with more than 2 dimensions')
                    return
                end
            end
            
            % set default parameters
            defaultX  = 1:size(signal,1);
            defaultY  = flip(1:size(signal,2));
            defaultSF = calc_scale_factor(signal);
            defaultC  = 'k';
            defaultLW = 1;
            defaultLS = '-';
            defaultZ = [];
            
            % parse inputs
            p = inputParser;
            addRequired(p,'ax1')
            addRequired(p,'signal');
            addParameter(p,'xvec',defaultX);
            addParameter(p,'yvec',defaultY);
            addParameter(p,'scalefactor',defaultSF)
            addParameter(p,'color',defaultC)
            addParameter(p,'linewidth',defaultLW)
            addParameter(p,'linestyle',defaultLS)
            addParameter(p,'zval',defaultZ)
            
            parse(p,ax1,signal,varargin{:});
            
            % set parameters
            xvec = p.Results.xvec;
            yvec = p.Results.yvec;
            scalefactor = p.Results.scalefactor;
            col = p.Results.color;
            linewidth = p.Results.linewidth;
            linestyle = p.Results.linestyle;
            zval = p.Results.zval;
            
            if numel(zval) == 1
                zval = zval*ones(size(signal,1),1);
            end
            
            % create figure if necessary
            if isempty(ax1)
                figure
                ax1 = axes('next','add');
            end
            set(ax1,'next','add')
            
            
            % plot
            nchans = size(signal,2);
            
            if isnumeric(col)
                if size(col,1) == 1
                    col = repmat(col,nchans,1);
                end
            else
                C    = cell(nchans,1);
                C(:) = {col};
                col  = C;
            end
            
            for j = 1:nchans
                chan_to_plot = signal(:,j) - nanmedian(signal(:,j));
                col4line     = col(j,:);
                if iscell(col4line)
                    col4line = col4line{1};
                end
                if isempty(zval)
                    plot(ax1,xvec,(chan_to_plot*scalefactor)+yvec(j),'color',col4line, 'Linewidth',linewidth,'Linestyle',linestyle);
                else
                    plot3(ax1,xvec,(chan_to_plot*scalefactor)+yvec(j),zval,'color',col4line, 'Linewidth',linewidth,'Linestyle',linestyle);
                end
            end
            
        end
        
        function drawNeuropixelSchematic(ax1, probe_type, varargin) % varargs = elec,shank,col,alpha
            
            if isempty(ax1)
                figure;
                ax1 = axes;
            end
            
            % Set up parameters
            [xvals, yvals] = get_x_and_y_vals(probe_type);
            
            refXvals = xvals(end-1:end); % reference sites for all probe types are currently even numbers (when 1-indexed), so on the right side of contacts
            internalRefs = NPadmin.getInteralRefElectrodes(probe_type);
            switch probe_type
                case 24
                    
                    nrows = 1280 ; % n electrodes = 1280; each electrode y=2 so nrows = n electrodes
                    totalShanks = 3; % index from 0
                    % 32 rows and in each block and 12 blocks per bank (last bank has only
                    % 4 blocks)
                    banklimits = [16*2*12 16*2*12*2 16*2*12*3];
                    
                case 21
                    nrows = 1280 ; % n electrodes = 1280; each electrode y=2 so nrows = n electrodes
                    totalShanks = 0; % index from 0
                    % 16 rows (2 electrodes per row) and in each block and 12 blocks per bank (last bank has only
                    % 4 blocks)
                    banklimits = [16*2*12 16*2*12*2 16*2*12*3];
                    
                case 0
                    nrows = 960 /2; % n electrodes = 960; but drawing 4 sites at a time rather than 2
                    totalShanks = 0; % index from 0
                    banklimits = [384 384*2];
                    
                    internalRefs = internalRefs -1; % not sure why but shift is needed to plot in correct place
            end
            
            % Draws the shanks
            % plot first four contacts then duplicate up whole probe
            % Matrix  goes clockwise then contacts left to right
            x = []; y = [];
            for i = 1:4:numel(xvals)
                x = [x xvals(i) xvals(i) xvals(i+1) xvals(i+1) xvals(i) NaN,...
                    xvals(i+2) xvals(i+2) xvals(i+3) xvals(i+3) xvals(i+2) NaN];
            end
            for i = 1:2:numel(yvals)
                y = [y yvals(i) yvals(i+1) yvals(i+1) yvals(i) yvals(i) NaN,...
                    yvals(i) yvals(i+1) yvals(i+1) yvals(i) yvals(i) NaN];
            end
            
            % Repeat for full length of probe
            if probe_type == 0
                ymultiplier = 0:4:(nrows*2)-1;
            else
                ymultiplier = 0:2:nrows-1;
            end
            ymultiplier = repelem(ymultiplier,numel(x));
            y = repmat(y,1,nrows/2);
            x = repmat(x,1,nrows/2);
            % Plot four times for each shank
            for i = 0:totalShanks
                % Plot electrodes
                hold(ax1,"on")
                plot(ax1,x+i,y+ymultiplier,'k','LineWidth',0.25');
                
                % Plot the bank limits
                for ii = 1:numel(banklimits)
                    plot(ax1,[min(xvals)-0.1 max(xvals)+0.1]+i,[banklimits(ii)  banklimits(ii)],'b','LineWidth',1.25)
                end
                
%                 % Plot internal references
%                 plot_square = @(row) plot(ax1,[refXvals(1) refXvals(1) refXvals(1) refXvals(2) refXvals(2) refXvals(2) refXvals(2) refXvals(1)]+i,...
%                     [row-1 row+1 row+1 row+1 row+1 row-1 row-1 row-1],'c','linewidth',1.5);
%                 for ii = 1:numel(internalRefs)
%                     plot_square(internalRefs(ii));
%                 end
                
                % Plot tip reference
                plot(ax1,[min(xvals) min(xvals)+((max(xvals)-min(xvals))/2)]+i,[0 -9],'c','LineWidth',1.5)
                plot(ax1,[min(xvals)+((max(xvals)-min(xvals))/2) max(xvals)]+i,[-9 0],'c','LineWidth',1.5)
                plot(ax1,[min(xvals) max(xvals)]+i,[0 0],'c','LineWidth',1.5)
                
            end
            
            ax1.XLabel.String = 'Shank'; ax1.YLabel.String = 'Electrode #';
            
            if nargin>2 % includes electrodes to shade
                elec  = varargin{1};
                shank = varargin{2};
                col   = varargin{3};
                alph  = varargin{4};
                
                for n=1:numel(elec)
                    NPplot.colourElectrodeSite(ax1,probe_type,shank(n),elec(n),col,alph)
                end
            end
            
        end
        function colourElectrodeSite(ax1,probe_type,shank,elec,col,alph)
            [xvals, ~] = get_x_and_y_vals(probe_type);
            if probe_type == 0
                % Need to have a slightly different system as there is a
                % checkerboard
                %                 [5] [6]
                %                   [3] [4]
                %                 [1] [2]
                if ismember(elec,1:4:2000) == 1 % Left column
                    x = [xvals(1)  xvals(1) xvals(2) xvals(2)];
                elseif ismember(elec,2:4:2000) == 1 % Right column
                    x = [xvals(3)  xvals(3) xvals(4) xvals(4)];
                elseif ismember(elec,3:4:2000) == 1 % staggered left column
                    x = [xvals(5)  xvals(5) xvals(6) xvals(6)];
                elseif ismember(elec,4:4:2000) == 1 % staggered left column
                    x = [xvals(7)  xvals(7) xvals(8) xvals(8)];
                end
            else
                if rem(elec,2) == 1 %Is odd then its on left side of shank
                    x = [xvals(1)  xvals(1) xvals(2) xvals(2)];
                else % Is even then on right side of shank
                    x = [xvals(3)  xvals(3) xvals(4) xvals(4)];
                end
            end
            x = x + shank-1; % Make adjustment for the shank
            row = elec;
            if ~rem(row,2)%is even
                row = row-1;
            end
            y = [row-0.9 row+0.9 row+0.9 row-0.9]; % each electrode covers 2 y units
            patch(ax1,x,y,col,'EdgeColor','none','FaceAlpha',alph);
        end
        function highlightElectrodeSite(ax1,probe_type,shank,elec,col)
            if isempty(elec)
                return
            end
             [xvals, ~] = get_x_and_y_vals(probe_type);
            if probe_type == 0
                % Need to have a slightly different system as there is a
                % checkerboard
                %                 [5] [6]
                %                   [3] [4]
                %                 [1] [2]
                if ismember(elec,1:4:2000) == 1 % Left column
                    x = [xvals(1)  xvals(1) xvals(2) xvals(2)];
                elseif ismember(elec,2:4:2000) == 1 % Right column
                    x = [xvals(3)  xvals(3) xvals(4) xvals(4)];
                elseif ismember(elec,3:4:2000) == 1 % staggered left column
                    x = [xvals(5)  xvals(5) xvals(6) xvals(6)];
                elseif ismember(elec,4:4:2000) == 1 % staggered left column
                    x = [xvals(7)  xvals(7) xvals(8) xvals(8)];
                end
            else
                if rem(elec,2) == 1 %Is odd then its on left side of shank
                    x = [xvals(1)  xvals(1) xvals(2) xvals(2)];
                else % Is even then on right side of shank
                    x = [xvals(3)  xvals(3) xvals(4) xvals(4)];
                end
            end
            x = x + shank-1; % Make adjustment for the shank
            row = elec;
            if ~rem(row,2)%is even
                row = row-1;
            end
            y = [row-0.9 row+0.9 row+0.9 row-0.9]; % each electrode covers 2 y units
            plot(ax1,[xvals(1) xvals(1) xvals(1) xvals(2) xvals(2) xvals(2) xvals(2) xvals(1)],...
                     [row-1 row+1 row+1 row+1 row+1 row-1 row-1 row-1],col,'linewidth',1.5);
              
        end
    end
end

function scalefactor = calc_scale_factor(signal) % subfunc for plotTraceAcrossChans
%posidx    = signal > 0;
% chans_neg = signal;
% chans_neg(posidx)=NaN;
% med_neg = nanmedian(chans_neg,1);
% chans_pos = signal;
% chans_pos(~posidx)=NaN;
% med_pos = nanmedian(chans_pos,1);
%meds = med_pos-med_neg;
meds = nanmedian(abs(signal));
max_diff = max(meds);
scalefactor = 1/max_diff;
end

function [xvals, yvals] = get_x_and_y_vals(probe_type)
if probe_type == 0
    xvals = [0.15 0.35 0.45 0.65 0 0.2 0.3 0.5];
    yvals = [0.1 1.9 2.1 3.9];
else
    xvals = [0 0.2 0.3 0.5];
    yvals = [0.1 1.9];
end
end
