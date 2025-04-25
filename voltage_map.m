clc
clearvars

% Deactivate warning for performance reasons
warning('off', 'signal:findpeaks:largeMinPeakHeight')

figure_dir = "figures";
data_dir = "processed_data";
patient_dirs = dir(data_dir);

n_patients = sum(~isnan(cellfun(@str2double, {patient_dirs.name})));
data = table('Size', [n_patients, 3], 'VariableTypes', ["string" "double" "double"], ...
    'VariableNames', ["Patient ID" "Entropy" "LASE"]);

% For each patient
for ipat = 1:length(patient_dirs)
    patient_dir = patient_dirs(ipat);
    
    % Debugging
    if patient_dir.name ~= "111"
        continue
    end
    
    % Skip the '.', '..' and '.DS_Store' directories
    if contains(patient_dir.name, '.')
        continue
    end

    % Load the data
    path_to_data = strjoin([data_dir patient_dir.name], '/');
    load(strjoin([path_to_data "LA_info.mat"], '/'))
    load(strjoin([path_to_data "LA_mesh.mat"], '/'))
    load(strjoin([path_to_data "LA_points_data.mat"], '/'))

    fs = double(fs);
   
    % Skip the patients with CS pacing
    if contains(map_name, "PACE")
        continue
    end
    
    disp("Processing patient " + patient_ID)

    columns = cellstr(columns);
    electrodes = cellstr(electrodes);
    signals = double(signals);
    triangles = double(triangles) + 1;
    
    % Check whether is AF or no
    AF = false;
    if contains(map_name, "FA") || contains(map_name, "AF")
        AF = true;
    end
    if patient_ID == "100"
        AF = true;
    end

    ecg_channel = 'V5';
    ecg_signal_index = find(strcmp(columns, ecg_channel));

    % Create empty matrices to store the PTP value and relative
    % coordinates
    coordinates = [];
    voltages = [];

    % Create the filters
    [b1, a1] = butter(3, (2*[30, 400])/fs);
    [b2, a2] = butter(3, (40)/fs);

    % For each point export
    for pp = 1:length(points_IDs)

        point_ID = points_IDs(pp);

        ecg_signal = signals(pp, :, ecg_signal_index);
        
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
        for ee = 1:length(electrodes)
            electrode_index = find(strcmp(electrodes{ee}, columns));
            
            egm_signal = signals(pp,:,electrode_index)';
            
            % Extract the windows from the signal
            egm_windows = arrayfun(@(a,b) egm_signal(a:b), ...
                egm_windows_bounds(:,1), egm_windows_bounds(:,2), 'Unif', 0);
            
            % Compute the peak to peak measure of each window
            egm_ptp = cellfun(@(x) peak2peak(x), egm_windows);

            % Get the coordinates of the electrode
            electrode_coordinates = squeeze(positions(pp, ee, :, :))';
            % Remove trailing 0s
            electrode_coordinates = electrode_coordinates( ...
                1:find(sum(electrode_coordinates, 2), 1, "last"),:);
            % Linearly interpolate coordinates to the length of the signals
            electrode_coordinates_timestamps = linspace( ...
                1, size(signals,2), length(electrode_coordinates));
            electrode_coordinates = interp1(electrode_coordinates_timestamps, ...
                electrode_coordinates, 1:length(signals), "linear");
            % Window the coordinates
            electrode_coordinates_windows = arrayfun(@(a,b) electrode_coordinates(a:b,:), ...
                egm_windows_bounds(:,1), egm_windows_bounds(:,2), "Unif", 0);
            % Compute the mean of the coordinates within the windows
            mean_electorde_coordinates = cellfun((@(x) mean(x,1)), ...
                electrode_coordinates_windows, "Unif", 0);
            mean_electorde_coordinates = cell2mat(mean_electorde_coordinates);

            coordinates = [coordinates; mean_electorde_coordinates]; %#ok<AGROW>
            voltages = [voltages; egm_ptp]; %#ok<AGROW>
        end
    end

    status = mkdir(strjoin([figure_dir, patient_ID], '/'));

    % % Sphere and cylinder computation on mesh
    % is_resampled = false;
    % vertex_voltage_map = vertex_voltage_mapping(vertices, triangles, voltages, coordinates, is_resampled);
    % final_voltage_map = max(vertex_voltage_map, [], 2);
    % 
    % % Plot mesh
    % figure()
    % title("Before resampling")
    % 
    % trisurf(triangles, vertices(:,1), vertices(:,2), vertices(:,3), ...
    %     final_voltage_map, 'edgecolor', 'none', 'facecolor', 'interp');
    % colormap(flipud(turbo));
    % colormap([0.5, 0.5, 0.5; flipud(turbo)]); % Append gray to the colormap
    % colorbar;
    % clb = colorbar;
    % clb.Limits = [0, 4];
    % material dull
    % cameraLight;
    % if AF, clim([0.05 0.24]); else, clim([0.05, 0.5]); end
    % savefig(strjoin([figure_dir, patient_ID, "voltage_map.fig"], '/'))
    % close
    % 
    % [f, entropy] = entropy_calculation(final_voltage_map);
    % savefig(f, strjoin([figure_dir, patient_ID, "entropy.fig"], '/'))
    % close

    % Write mesh to disk for meshtool
    vtkwrite(strjoin([data_dir, patient_ID, 'LA_mesh.vtk'], '/'), 'polydata', ...
         'triangle', vertices(:,1), vertices(:,2), vertices(:,3), triangles);

    % Resample mesh with meshtool
    command = sprintf("./meshtool/meshtool resample surfmesh -msh=%s -avrg=5 " + ...
        "-outmsh=%s -ofmt=vtk_polydata -surf_corr=0.8", ...
        strjoin([data_dir, patient_ID, 'LA_mesh.vtk'], '/'), ...
        strjoin([data_dir, patient_ID, 'LA_mesh_resampled.vtk'], '/'));
    [status, ~] = system(command);
    if status ~= 0
        error("Meshtool exit status is non-zero.")
    end

    [vertices_rsmp, triangles_rsmp] = read_vtk(strjoin([data_dir, patient_ID, 'LA_mesh_resampled.vtk'], '/'));
    vertices_rsmp = vertices_rsmp';
    triangles_rsmp = triangles_rsmp';

    % Sphere and cylinder computation on mesh
    is_resampled = true;
    vertex_voltage_map_rsmp = vertex_voltage_mapping(vertices_rsmp, triangles_rsmp, voltages, coordinates, is_resampled);
    final_voltage_map_rsmp = max(vertex_voltage_map_rsmp, [], 2);
    
    % Plot resampled mesh
    figure()
    title("After resampling")
    trisurf(triangles_rsmp, vertices_rsmp(:,1), vertices_rsmp(:,2), vertices_rsmp(:,3), ...
        final_voltage_map_rsmp, 'edgecolor', 'none', 'facecolor', 'interp');
    colormap(flipud(turbo));
    colormap([0.5, 0.5, 0.5; flipud(turbo)]); % Append gray to the colormap
    colorbar;
    clb = colorbar;
    clb.Limits = [0, 4];
    material dull
    cameraLight;
    if AF, clim([0.05 0.24]); else, clim([0.05, 0.5]); end
    savefig(strjoin([figure_dir, patient_ID, "voltage_map_rsmp.fig"], '/'))
    close
    
    [f_rsmp, lase] = entropy_calculation(final_voltage_map_rsmp);
    savefig(f_rsmp, strjoin([figure_dir, patient_ID, "entropy_rsmp.fig"], '/'))
    close

    data(ipat,:) = {patient_ID, 0, lase};
end
writetable(data, "entropy.csv")
