function [roiVoxelXYZ, roiCoreWorldXYZ] = f_xyz2roi_withCoreg(oriSpaceMRI, tarSpaceMRI, oriWorldXYZ, Radi)


%% main function
    % Co-registration the original space T1w to the target space T1w image
    cfg = [];
    cfg.method = 'spm';
    cfg.spmversion = 'spm12';
    cfg.coordsys = 'acpc';
    cfg.viewresult = 'yes';
    ori2tarMRI = ft_volumerealign(cfg, oriSpaceMRI, tarSpaceMRI);
    % Double check the co-registration
    if sum(sum(sum(ori2tarMRI.anatomy - oriSpaceMRI.anatomy))) ~= 0
        error("wrong co-registration"); 
    end

    % Convert original space world coordinate to ROI in target space
    oriAffine = oriSpaceMRI.transform;
    ori2tarAffine = ori2tarMRI.transform;
    for n = 1:size(oriWorldXYZ, 1)
        oriVoxelXYZ = f_world2voxel(oriWorldXYZ(n, :), oriAffine); % obtain the coordiante of each channel in voxcel space
        tarWorldXYZ = f_voxel2world(oriVoxelXYZ, ori2tarAffine); % obtain the coodinate of each channel in atlas reference space
        roiVoxelXYZ{n} = f_xyz2roi(tarWorldXYZ, tarSpaceMRI, Radi);
        roiCoreWorldXYZ(n, :) = tarWorldXYZ;
    end
end