function vertex_voltage_map = vertex_voltage_mapping(vertices, triangles, voltages, coordinates, is_resampled, thrs)

    % Keep only voltages less than given threshold
    voltage_upper_bound = thrs;
    coordinates = coordinates(voltages < voltage_upper_bound, :);
    voltages = voltages(voltages < voltage_upper_bound);
    
    r = 8;
    d = 2;

    TR = triangulation(triangles, vertices);
    vertices_normals = vertices + vertexNormal(TR);

    if is_resampled
        E = edges(TR);
        avg_edge_length = mean(sqrt(sum((TR.Points(E(:,1), :) - TR.Points(E(:,2), :)).^2, 2)));
        d = 2 * avg_edge_length;
    end
    
    % Sphere
    sphere_candidates = pdist2(vertices, coordinates, 'squaredeuclidean') <= r^2;
    
    % Cylinder
    cylinder_distance = @(u, v, n_v) norm(cross(u - v, u - n_v)); % / norm(n_v) == 1;
    cylinder_candidates_distances = ones(size(sphere_candidates)) + d;
    [vertices_idx, coordinates_idx] = find(sphere_candidates);
    for idx = 1:length(vertices_idx)
        ii = vertices_idx(idx);
        jj = coordinates_idx(idx);

        u = vertices(ii,:);
        v = coordinates(jj,:);
        n_v = vertices_normals(ii,:);

        cylinder_candidates_distances(ii,jj) = cylinder_distance(u, v, n_v);
    end
    final_candidates = cylinder_candidates_distances <= d;
    
    vertex_voltage_map = voltages' .* final_candidates;
end
