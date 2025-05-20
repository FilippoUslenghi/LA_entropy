clc
clearvars

data_dir = "processed_data";
patient_dirs = dir(data_dir);
out_dir = "results";

experiments = ["thrs_<15_no_filt",  "thrs_<15_filt"];
thrss = [15, 15];

verbose = false;
for iexp = 1:length(experiments)
    thrs = thrss(iexp);
    experiment = experiments(iexp);
    mkdir(strjoin([out_dir, experiment], '/'))
    
    disp([experiment, thrs, iexp])
    
    n_patients = sum(~isnan(cellfun(@str2double, {patient_dirs.name})));
    data = table('Size', [n_patients, 3], 'VariableTypes', ["string" "double" "double"], ...
        'VariableNames', ["Patient ID" "Entropy" "LASE"]);
    
    patients_rhythms = load("patients_rhythms.mat").patients_rhythms;

    % For each patient
    parfor ipat = 1:length(patient_dirs)
        % Deactivate warning for performance reasons
        warning('off', 'signal:findpeaks:largeMinPeakHeight')
        warning('off', 'MATLAB:MKDIR:DirectoryExists')

        patient_dir = patient_dirs(ipat);
        
        % Skip the '.', '..' and '.DS_Store' directories
        if contains(patient_dir.name, '.')
            continue
        end
        % Skip the patient 118 because it has too many points
        if patient_dir.name == "118"
            data(ipat,:) = {"118", 0, 0};
            continue
        end
    
        % Load the data
        path_to_data = strjoin([data_dir patient_dir.name], '/');
        INFO = load(strjoin([path_to_data "LA_info.mat"], '/'))
        MESH = load(strjoin([path_to_data "LA_mesh.mat"], '/'))
        POINTS = load(strjoin([path_to_data "LA_points_data.mat"], '/'))
    
        INFO.fs = double(INFO.fs);
        
        disp("Processing patient " + INFO.patient_ID + newline)
    
        POINTS.columns = cellstr(POINTS.columns);
        POINTS.electrodes = cellstr(POINTS.electrodes);
        POINTS.signals = double(POINTS.signals);
        MESH.triangles = double(MESH.triangles) + 1;
        
        % Check whether is AF or no
        AF = false;
        rhythm = string(patients_rhythms{strcmp(patients_rhythms(:,1), '100'), 2});
        if rhythm ~= "SR"
            AF = true;
        end
    
        ecg_channel = 'V5';
        ecg_signal_index = find(strcmp(POINTS.columns, ecg_channel));
    
        % Create empty matrices to store the PTP value and relative
        % coordinates
        coordinates = [];
        voltages = [];
    
        % Create the filters
        [b1, a1] = butter(3, (2*[30, 300])/INFO.fs);
        [b2, a2] = butter(3, (40)/INFO.fs);
    
        % For each point export
        for pp = 1:length(POINTS.points_IDs)
    
            point_ID = POINTS.points_IDs(pp);
    
            ecg_signal = POINTS.signals(pp, :, ecg_signal_index);
            
            % Apply first filter
            filtered_ecg_signal = filtfilt(b1, a1, ecg_signal);
            abs_filtered_ecg_signal = abs(filtered_ecg_signal);
            
            % Apply second filter
            ecg_activation = filtfilt(b2, a2, abs_filtered_ecg_signal);
    
            thr = 0.08;
            [~, ecg_peak_train] = findpeaks(ecg_activation, "MinPeakHeight", thr, ...
                "MinPeakDistance", 600);
    
            % Create the windows
            window_bounds = [-50; 100];
            ecg_windows = ecg_peak_train + window_bounds;
            ecg_windows = min(max(1,ecg_windows), 2500);
    
            % Exclude the ecg windows from the EGM signals
            egm_windows_bounds = [ecg_windows, [length(ecg_signal); 1]];
            egm_windows_bounds = reshape(circshift(egm_windows_bounds(:), 1), ...
                size(egm_windows_bounds))';
    
            % For each electrode
            for ee = 1:length(POINTS.electrodes)
                electrode_index = find(strcmp(POINTS.electrodes{ee}, POINTS.columns));
                
                egm_signal = POINTS.signals(pp,:,electrode_index)';
                if mod(iexp,2) == 0
                    egm_signal = filtfilt(b1, a1, egm_signal);
                end
    
                % Extract the windows from the signal
                egm_windows = arrayfun(@(a,b) egm_signal(a:b), ...
                    egm_windows_bounds(:,1), egm_windows_bounds(:,2), 'Unif', 0);
                
                % Compute the peak to peak measure of each window
                egm_ptp = cellfun(@(x) peak2peak(x), egm_windows);
    
                % Get the coordinates of the electrode
                electrode_coordinates = squeeze(POINTS.positions(pp, ee, :, :))';
                % Remove trailing 0s
                electrode_coordinates = electrode_coordinates( ...
                    1:find(sum(electrode_coordinates, 2), 1, "last"),:);
                % Linearly interpolate coordinates to the length of the signals
                electrode_coordinates_timestamps = linspace( ...
                    1, size(POINTS.signals,2), length(electrode_coordinates));
                electrode_coordinates = interp1(electrode_coordinates_timestamps, ...
                    electrode_coordinates, 1:length(POINTS.signals), "linear");
                % Window the coordinates
                electrode_coordinates_windows = arrayfun(@(a,b) electrode_coordinates(a:b,:), ...
                    egm_windows_bounds(:,1), egm_windows_bounds(:,2), "Unif", 0);
                % Compute the mean of the coordinates within the windows
                mean_electorde_coordinates = cellfun((@(x) mean(x,1)), ...
                    electrode_coordinates_windows, "Unif", 0);
                mean_electorde_coordinates = cell2mat(mean_electorde_coordinates);
    
                coordinates = [coordinates; mean_electorde_coordinates];
                voltages = [voltages; egm_ptp];
            end
        end
    
        % Sphere and cylinder computation on mesh
        is_resampled = false;
        vertex_voltage_map = vertex_voltage_mapping(MESH.vertices, MESH.triangles, voltages, coordinates, is_resampled, thrs);
        final_voltage_map = max(vertex_voltage_map, [], 2);
        [f, entropy] = entropy_calculation(final_voltage_map, verbose);
    
        % Write mesh to disk for meshtool
        vtkwrite(strjoin([data_dir, INFO.patient_ID, 'LA_mesh.vtk'], '/'), 'polydata', ...
             'triangle', MESH.vertices(:,1), MESH.vertices(:,2), MESH.vertices(:,3), MESH.triangles);
    
        % Resample mesh with meshtool
        command = sprintf("./meshtool/meshtool resample surfmesh -msh=%s -avrg=5 " + ...
            "-outmsh=%s -ofmt=vtk_polydata -surf_corr=0.8", ...
            strjoin([data_dir, INFO.patient_ID, 'LA_mesh.vtk'], '/'), ...
            strjoin([data_dir, INFO.patient_ID, 'LA_mesh_resampled.vtk'], '/'));
        [status, ~] = system(command);
        if status ~= 0
            error("Meshtool exit status is non-zero.")
        end
    
        [vertices_rsmp, triangles_rsmp] = read_vtk(strjoin([data_dir, INFO.patient_ID, 'LA_mesh_resampled.vtk'], '/'));
        vertices_rsmp = vertices_rsmp';
        triangles_rsmp = triangles_rsmp';
    
        % Sphere and cylinder computation on mesh
        is_resampled = true;
        vertex_voltage_map_rsmp = vertex_voltage_mapping(vertices_rsmp, triangles_rsmp, voltages, coordinates, is_resampled, thrs);
        final_voltage_map_rsmp = max(vertex_voltage_map_rsmp, [], 2);
        [f_rsmp, lase] = entropy_calculation(final_voltage_map_rsmp, verbose);
    
        data(ipat,:) = {INFO.patient_ID, 0, lase};
    end
    writetable(data, strjoin([out_dir, experiment, "lase.csv"], '/'))
end
