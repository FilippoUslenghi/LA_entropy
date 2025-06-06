% Script that downsamples the original meshes created by CARTO 3 system
% using meshtool
clc
clearvars

data_dir = "processed_data";
dir_struct = dir(data_dir);
num_vertices = [];
for i = 1:size(dir_struct, 1)
    patient_dir = dir_struct(i);
    patient_ID = patient_dir.name;
    
    % Skip directory such as '.', '..' and '.DS_Store'
    if contains(patient_ID, ".")
        continue
    end

    mesh_path = strjoin([data_dir, patient_ID, ], '/');

    command = sprintf("./meshtool/meshtool resample surfmesh -msh=%s -avrg=2 " + ...
    "-outmsh=%s -ofmt=vtk_polydata -surf_corr=0.8", ...
    strjoin([data_dir, patient_ID, 'LA_mesh.vtk'], '/'), ...
    strjoin([data_dir, patient_ID, 'LA_mesh_resampled.vtk'], '/'));
    [status, ~] = system(command);
    if status ~= 0
        error("Meshtool exit status is non-zero.")
    end

    [vertices, triangles] = read_vtk(strjoin([data_dir, ...
        patient_ID, 'LA_mesh_resampled.vtk'], '/'));
    num_vertices = [num_vertices; size(vertices,2)];
end
