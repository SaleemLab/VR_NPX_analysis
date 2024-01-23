%%facemap file handling
base_folder = '/mnt/rds01/ibn-vision/DATA/SUBJECTS/';
mouse = {'M23087','M23031','M23032','M23034','M23037','M23038','M23017','M23028','M23029','M23153'};
excluded_dirs = {'20231214', '20231214_OpenField', '20231215', '20231215_OpenField'};

for i = 1:length(mouse)
    mouse_id = mouse{i};
    date_dir_ephys = dir(fullfile(base_folder, mouse_id, 'stimuli'));
    date_dir_ephys = {date_dir_ephys([date_dir_ephys.isdir]).name};
    date_dir_ephys = date_dir_ephys(~ismember(date_dir_ephys, {'.', '..'}));
    
    for j = 1:length(date_dir_ephys)
        date_dir = date_dir_ephys{j};
        stimulus_folder = fullfile(base_folder, mouse_id, 'stimuli', date_dir);
        
        if ~ismember(date_dir, excluded_dirs)
            files = dir(fullfile(stimulus_folder, '*'));
            files = {files.name};
            
            avi_files = files(contains(files, '.avi') & ~contains(files, 'shuffle'));
            disp(avi_files);
            
            proc_files = strrep(avi_files, '.avi', '_proc.mat');
            disp(proc_files);
            
            initial_proc_file = files(find(contains(files, '_proc.mat'), 1));
            
            for k = 1:length(proc_files)
                file = proc_files{k};
                if ~strcmp(file, initial_proc_file{1})
                    copyfile(fullfile(stimulus_folder, initial_proc_file{1}), fullfile(stimulus_folder, file));
                end
            end
            
            for k = 1:length(proc_files)
                proc_file = proc_files{k};
                avi_file = avi_files{k};
                if ~strcmp(proc_file, initial_proc_file{1})
                    proc = load(fullfile(stimulus_folder, proc_file));
                    filenames_tmp = fullfile(stimulus_folder, avi_file);
                    proc.filenames = reshape(filenames_tmp,[1,1,size(filenames_tmp,2)]);
                    save(fullfile(stimulus_folder, proc_file), '-struct', 'proc');
                end
            end
        end
    end
end