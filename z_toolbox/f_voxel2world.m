function worldxyz = f_voxel2world(voxelxyz, Affine)
    % convert the coordinate of voxel space to the reference
    worldxyz = Affine * [voxelxyz, ones(length(voxelxyz(:,1)),1)]' ;
    worldxyz = worldxyz';
    worldxyz(:, 4) = [];
end