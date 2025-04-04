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

fs = 1000;
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
   
    % Skip the patients with CS pacing
    if contains(map_name, "PACE")
        continue
    end
    
    disp("Processing patient " + patient_ID)

    columns = cellstr(columns);
    electrodes = cellstr(electrodes);
    signals = double(signals);
    
    % Select reference channel
    AF = true;
    if ~contains(map_name, "FA")
        AF = false;
        reference_channel = 'CS1-CS2';
    end
    if patient_ID == "100"
        AF = true;
        reference_channel = 'V5';
    end

    disp("Reference channel is " + reference_channel)
    reference_signal_index = find(strcmp(columns, reference_channel));

    % For each point export
    for pp = 1:length(points_IDs)
        point_ID = points_IDs(pp);

        reference_signal = signals(pp, :, reference_signal_index);
        
        % First filter
        [b, a] = butter(3, (2*[30, 400])/fs);
        filtered_reference_signal = filtfilt(b, a, reference_signal);
        % Second filter
        [b, a] = butter(3, (40)/fs);
        abs_filtered_reference_signal = abs(filtered_reference_signal);

        reference_activation = filtfilt(b, a, abs_filtered_reference_signal);

        thr = 0.08;
        [pks, peak_train] = findpeaks(reference_activation, 'MinPeakHeight', ...
            thr,'MinPeakDistance',600);
        
        window_bounds = [-150 30];
        windows = peak_train' + window_bounds;

        if AF
            % Shift on the left the windows by a certain amount
            windows = windows - 100;
        end
        windows = max(windows,1);

        % For each electrode
        for ee = 1:length(electrodes)
            electrode_index = find(strcmp(electrodes{ee}, columns));
            
            egm_signal = signals(pp,:,electrode_index);
            
            % Extract the windows from the signal
            egm_windows = arrayfun((@(a,b) egm_signal(a:b)), ...
                windows(:,1), windows(:,2), 'Unif', 0);
            
            % Compute the peak to peak measure of each window
            egm_ptp = cellfun(@(x) peak2peak(x), egm_windows);

            % TODO: find the relative positions of the EGM activations
            % (tip: compute the mean of the coordinates in each time window)
            

            break
        end

        % TODO: create a table that associates each positions
        % to its relative EGM activation


        break
    end
    
    break
end












