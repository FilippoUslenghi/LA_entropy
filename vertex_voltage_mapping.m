function vertex_voltage_map = vertex_voltage_mapping(vertices, triangles, voltages, coordinates, is_resampled)

    r = 8;
    d = 2;

    % if is_resampled
    %     % Get neighboring vertices    
    %     neighbors = zeros(size(Vertices,1), size(Vertices,1));
    %     for ii = 1:size(Vertices,1)
    %         [row,col] = find(Triangles==ii);
    %         for jj = 1:length(row)
    %             if col(jj)==1
    %                 neighbors(ii, Triangles(row(jj),2)) = 1;
    %                 neighbors(ii, Triangles(row(jj),3)) = 1;
    %             end
    %             if col(jj)==2
    %                 neighbors(ii, Triangles(row(jj),1)) = 1;
    %                 neighbors(ii, Triangles(row(jj),3)) = 1;
    %             end
    %             if col(jj)==3
    %                 neighbors(ii, Triangles(row(jj),1)) = 1;
    %                 neighbors(ii, Triangles(row(jj),2)) = 1;
    %             end
    %         end
    %     end
    % 
    %     % Distance between vertices.
    %     distances = zeros(size(Vertices,1), size(Vertices,1));
    %     for ii=1:size(Vertices,1)
    %         for jj=1:size(Vertices,1)
    %             if neighbors(ii, jj) == 1
    %                 distances(ii,jj) =  pdist([Vertices(ii,:); Vertices(jj,:)], 'euclidean'); %distanze solo dei vicini
    %             end
    %         end
    %     end
    %     d = 2 * mean(distances(find(distances)));
    
    vertices_normals = vertices + vertexNormal(triangulation(triangles, vertices));
    
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

