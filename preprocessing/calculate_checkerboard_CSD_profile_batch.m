function calculate_checkerboard_CSD_profile_batch(experiment_info,Stimulus_type)
% calculate_checkerboard_CSD_profile_batch - 
%
%   calculate_checkerboard_CSD_profile_batch(experiment_info,Stimulus_type)
%
%   This function does [brief description of what the function does].
%
%   INPUTS:
%   - input1: [datatype] Description of the first input parameter.
%   - input2: [datatype] Description of the second input parameter.
%
%   OUTPUTS:
%   - output1: [datatype] Description of the first output parameter.
%   - output2: [datatype] Description of the second output parameter.
%
%   NOTES:
%   - Additional notes or important information about the function.
%
%   EXAMPLE:
%   % Example usage of the function
%   [result1, result2] = myFunction(input1_value, input2_value);
%
%   SEE ALSO:
%   - Related functions or scripts
%
%   AUTHOR:
%   Your Name
%
%   DATE:
%   Current date

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    
    if isempty(stimulus_name)
        fprintf('%s %s missing checkerboard',session_info.probe(1).SUBJECT,session_info.probe(1).SESSION);
        continue
    end

    for n = 1:length(session_info) % Just in case if there are multiple sessions from the same stimulus (should retire in the future)
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options= session_info(n).probe(nprobe);
%             options.ROOTPATH = ROOTPATH;
            options.importMode = 'LF';
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            [lfpAvg.nprobe(options.probe_no),csd.nprobe(options.probe_no),power,best_channels] = checkerboard_CSD_profile(options);
            save_all_figures(options.ANALYSIS_DATAPATH,[]);
            close
            
            save(fullfile(options.ANALYSIS_DATAPATH,"checkerboard_CSD.mat"),'lfpAvg','csd');

            column = 1;
            [LF_FILE imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
            [best_channels_updated{options.probe_no}] = update_best_channels(options,sorted_config);
            
            % If updating all
            best_channels{options.probe_no}.first_in_brain_channel = best_channels_updated{options.probe_no}(1);
            best_channels{options.probe_no}.L4_channel = best_channels_updated{options.probe_no}(2);
%             best_channels{nprobe}.L4_channel = [];
            best_channels{options.probe_no}.L5_channel = best_channels_updated{options.probe_no}(3);
            best_channels{options.probe_no}.CA1_channel = best_channels_updated{options.probe_no}(4);
            
            close
            
            save(fullfile(erase(options.ANALYSIS_DATAPATH,['\','Checkerboard']),"best_channels.mat"),'best_channels')
            % Replot based on updated channels
            plot_perievent_CSD_LFP(lfpAvg.nprobe(options.probe_no),csd.nprobe(options.probe_no),power{options.probe_no},chan_config,sorted_config,best_channels{options.probe_no},options)
           save_all_figures(options.ANALYSIS_DATAPATH,[]);
        end
    end
end

