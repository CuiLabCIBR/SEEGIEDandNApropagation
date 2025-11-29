function voxelxyz = f_world2voxel(worldxyz, Affine)
    % convert the coordinate of reference space to voxel space 
    voxelxyz = inv(Affine)*[worldxyz, ones(length(worldxyz(:, 1)),1)]';
    voxelxyz = voxelxyz';
    voxelxyz(:, 4) = [];
end