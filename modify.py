import sys
import re

with open('main_extract_LFP_sleep.m', 'r') as f:
    code = f.read()

# Replace initialization block
init_pattern = r'for nprobe = 1:2\n\s*periripple_LFP_info_V1\(nprobe\)\.SO_amplitude = \[\];.*?end'
init_replacement = '''for nprobe = 1:2
    % periripple_LFP_info_V1(nprobe).SO_amplitude = [];
    % periripple_LFP_info_V1(nprobe).SO_amplitude{2} = [];
    
    % Initialize UP and DOWN onset arrays
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_amplitude{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_amplitude{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_amplitude{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_amplitude{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_phase{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_phase{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_phase{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_phase{2} = [];

    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_amplitude{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_amplitude{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_amplitude{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_amplitude{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_phase{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_phase{2} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_phase{1} = [];
    periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_phase{2} = [];

    periUP_DOWN_LFP_info_V1(nprobe).tvec = -1:0.02:1;
end'''
code = re.sub(init_pattern, init_replacement, code, flags=re.DOTALL)


# Replace the main loop block
loop_pattern = r'for nprobe = 1:2\n\s*event_sample_idx = round\(interp1\(lfp_V1\.timestamps.*?end\n\s*end'
loop_replacement = '''for nprobe = 1:2
        % --- Commented out ripple extraction ---
        % event_sample_idx = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), ripples(nprobe).onset, 'nearest', 'extrap'));
        % for nevent = 1:length(event_sample_idx)
        %     ripple_windows = event_sample_idx(nevent) + win_full_samp; % ripple windows
        %     ... (ripple and spindle extraction commented out) ...
        % end

        % --- Extract for UP onset ---
        UP_onset_times = slow_waves(nprobe).ints.UP(:,1);
        event_sample_idx_UP = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), UP_onset_times, 'nearest', 'extrap'));
        
        for nevent = 1:length(event_sample_idx_UP)
            event_windows = event_sample_idx_UP(nevent) + win_full_samp; 
            if sum(event_windows > length(lfp_V1.SO_amplitude_LFP))>0
                event_windows(sum(event_windows <= length(lfp_V1.SO_amplitude_LFP)):end) = nan;
            elseif sum(event_windows < 1)>0
                event_windows(event_windows<1) = nan;
            end
            valid_bins = ~isnan(event_windows);

            for mprobe = 1:2
                % SO amplitude
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_amplitude{mprobe} output];
                
                % Spindle amplitude
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_amplitude{mprobe} output];
                
                % SO phase
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_onset.SO_phase{mprobe} output];
                
                % Spindle phase
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_onset.spindle_phase{mprobe} output];
            end
        end

        % --- Extract for DOWN onset (UP offset) ---
        DOWN_onset_times = slow_waves(nprobe).ints.DOWN(:,1);
        event_sample_idx_DOWN = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), DOWN_onset_times, 'nearest', 'extrap'));
        
        for nevent = 1:length(event_sample_idx_DOWN)
            event_windows = event_sample_idx_DOWN(nevent) + win_full_samp; 
            if sum(event_windows > length(lfp_V1.SO_amplitude_LFP))>0
                event_windows(sum(event_windows <= length(lfp_V1.SO_amplitude_LFP)):end) = nan;
            elseif sum(event_windows < 1)>0
                event_windows(event_windows<1) = nan;
            end
            valid_bins = ~isnan(event_windows);

            for mprobe = 1:2
                % SO amplitude
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_amplitude{mprobe} output];
                
                % Spindle amplitude
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_amplitude{mprobe} output];
                
                % SO phase
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.SO_phase{mprobe} output];
                
                % Spindle phase
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_onset.spindle_phase{mprobe} output];
            end
        end
    end'''
code = re.sub(loop_pattern, loop_replacement, code, flags=re.DOTALL)


code = code.replace("save(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'),'periripple_LFP_info_V1_best_SO','-v7.3');", "% save(fullfile(analysis_folder,'periripple_LFP_info_V1_best_SO.mat'),'periripple_LFP_info_V1_best_SO','-v7.3');")
code = code.replace("save(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'),'periripple_LFP_info_V1','-v7.3');", "% save(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'),'periripple_LFP_info_V1','-v7.3');\n    save(fullfile(analysis_folder,'periUP_DOWN_LFP_info_V1.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');")
code = code.replace("save(fullfile(analysis_folder,'periripple_LFP_info_V1_PRE.mat'),'periripple_LFP_info_V1','-v7.3');", "% save(fullfile(analysis_folder,'periripple_LFP_info_V1_PRE.mat'),'periripple_LFP_info_V1','-v7.3');\n    save(fullfile(analysis_folder,'periUP_DOWN_LFP_info_V1_PRE.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');")
code = code.replace("save(fullfile(analysis_folder,'periripple_LFP_info_V1_RUN.mat'),'periripple_LFP_info_V1','-v7.3');", "% save(fullfile(analysis_folder,'periripple_LFP_info_V1_RUN.mat'),'periripple_LFP_info_V1','-v7.3');\n    save(fullfile(analysis_folder,'periUP_DOWN_LFP_info_V1_RUN.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');")


with open('main_extract_LFP_sleep.m', 'w') as f:
    f.write(code)

print('Success')
