% Scarta i pazienti con pacing del CS
% Considera sempre il CS, e.g. CS1-CS2 come reference
% Filtra i segnali col codice di Francesco per trovare le attivazioni
% Centra una finestra nel massimo [dell'attivazione del CS filtrato, escludendo il QRS' ...
% e.g.[x-120, x+30].
% Applica questa finestra a tutti i segnali degli EGM (cioè considera tutti gli elettrodi)
% Calcola il PTP di ogni finestra di EGM e associa ad ogni PTP
% la posizione (x,y,z) dell'elettrodo in quell'istante di tempo
% Crea una tabella finale che associ i valori di PTP alle sue coordinate 3D

% Per ogni punto della mesh, applica il metodo della sfera e del cilindro
% Considerando i punti della tabella.
clc
clearvars

% Deactivate warning for performance reasons
warning('off', 'signal:findpeaks:largeMinPeakHeight')

data_dir = "processed_data";
patient_dirs = dir(data_dir);

% For each patient
for ipat = 1:length(patient_dirs)
    patient_dir = patient_dirs(ipat);

    % Skip the '.' and '..' directories
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

    reference_channel = 'CS1-CS2';
    reference_signal_index = find(strcmp(columns, reference_channel));

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

        reference_signal = signals(pp, :, reference_signal_index);
        ecg_signal = signals(pp, :, ecg_signal_index);
        
        % Apply first filter
        filtered_reference_signal = filtfilt(b1, a1, reference_signal);
        filtered_ecg_signal = filtfilt(b1, a1, ecg_signal);
        
        abs_filtered_reference_signal = abs(filtered_reference_signal);
        abs_filtered_ecg_signal = abs(filtered_ecg_signal);
        
        % Apply second filter
        reference_activation = filtfilt(b2, a2, abs_filtered_reference_signal);
        ecg_activation = filtfilt(b2, a2, abs_filtered_ecg_signal);

        thr = 0.08;
        [~, reference_peak_train] = findpeaks(reference_activation, ...
            'MinPeakHeight', thr,'MinPeakDistance',600);
        [~, ecg_peak_train] = findpeaks(ecg_activation, "MinPeakHeight", thr, ...
            "MinPeakDistance", 600);

        % Create the windows
        window_bounds = [-100 100];
        reference_windows = reference_peak_train' + window_bounds;
        ecg_windows = ecg_peak_train' + window_bounds;
        
        % Exclude the ecg windows from the reference windows, works if
        % windows have same length
        for rr = 1:size(reference_windows, 1)
            reference_window = reference_windows(rr, :);
            for ee = 1:size(ecg_windows, 1)
                ecg_window = ecg_windows(ee, :);
                
                if min([reference_window ecg_window]) == min(reference_window)
                    reference_window(2) = min([reference_window(2) ecg_window]);
                else
                   reference_window(1) = max([reference_window(1) ecg_window]);
                end
            end
            % If reference window has length 0, mark to discard
            if min(reference_window) == max(reference_window)
                reference_window = [-1 -1];
            end
            reference_windows(rr, :) = reference_window;
        end

        % Discard previously marked windows
        reference_windows = reference_windows(any(reference_windows>1, 2), :);

        % Crop windows that exceed the boundaries
        reference_windows = max(reference_windows,1);
        reference_windows = min(reference_windows, 2500);

        % For each electrode
        for ee = 1:length(electrodes)
            electrode_index = find(strcmp(electrodes{ee}, columns));
            
            egm_signal = signals(pp,:,electrode_index)';
            
            % Extract the windows from the signal
            egm_windows = arrayfun(@(a,b) egm_signal(a:b), ...
                reference_windows(:,1), reference_windows(:,2), 'Unif', 0);
            
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
                reference_windows(:,1), reference_windows(:,2), "Unif", 0);
            % Compute the mean of the coordinates within the windows
            mean_electorde_coordinates = cellfun((@(x) mean(x,1)), ...
                electrode_coordinates_windows, "Unif", 0);
            mean_electorde_coordinates = cell2mat(mean_electorde_coordinates);

            coordinates = [coordinates; mean_electorde_coordinates]; %#ok<AGROW>
            voltages = [voltages; egm_ptp]; %#ok<AGROW>
        end
    end

    % Write mesh to disk for MeshTool
    % vtkwrite(strcat(data_dir, '/', patient_ID, '/', 'LA_mesh.vtk'), 'polydata', ...
    %      'triangle', vertices(:,1), vertices(:,2), vertices(:,3), triangles);

    % Resample mesh with MeshTool
    % command = ['./MeshTool resample surfmesh -msh=LA_mesh.vtk -avrg=5 ' ...
    %     '-outmsh=resampled_LA_mesh.vtk -ofmt=vtk_polydata -surf_corr=0.8'];
    % status = system(command);
    % if status ~= 0
    %     error("MeshTool exit status is non-zero.")
    % end

    % Sphere and cylinder computation on mesh
    is_resampled = false;
    vertex_voltage_map = vertex_voltage_mapping(vertices, triangles, voltages, coordinates, is_resampled);
    final_voltage_map = max(vertex_voltage_map, [], 2);
    break
end
