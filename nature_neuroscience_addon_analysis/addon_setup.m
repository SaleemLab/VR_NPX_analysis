function rev_dir = addon_setup()
% Path and data-folder setup for the Nature Neuroscience revision analyses,
% called at the top of every main_*.m.
% Adds this folder and all its subfolders to the MATLAB path and returns the
% revision analysis folder:
%   <D: or P:>\corticohippocampal_replay\nature_neuroscience_revision_analysis
% The rest of VR_NPX_analysis must also be on the path (e.g.
% subject_session_stimuli_mapping, merge_bilateral_ripple_events).
    root = fileparts(mfilename('fullpath'));
    addpath(genpath(root));
    if isempty(which('merge_bilateral_ripple_events'))
        warning('addon_setup:repo', ['VR_NPX_analysis is not on the MATLAB path; add it with ' ...
            'addpath(genpath(''%s''))'], fileparts(root));
    end
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    else
        error('addon_setup:data', 'corticohippocampal_replay folder not found on D: or P:');
    end
    rev_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis');
end
