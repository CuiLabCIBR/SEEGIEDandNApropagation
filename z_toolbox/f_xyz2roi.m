function roiVoxelXYZ = f_xyz2roi(worldxyz, Image, Radi)

%% main function
    Affine = Image.transform;
    imageSize = Image.dim;
    voxelxyz = f_world2voxel(worldxyz, Affine);

    ResoX = Image.hdr.xsize; RangeX = ceil(Radi / ResoX)+1;
    ResoY = Image.hdr.ysize; RangeY = ceil(Radi / ResoY)+1;
    ResoZ = Image.hdr.zsize; RangeZ = ceil(Radi / ResoZ)+1;

    roiSize = 0; roiVoxelXYZ = [];
    for dx = -RangeX : RangeX
        for dy = -RangeY : RangeY
            for dz = -RangeZ : RangeZ
                X = voxelxyz(1) + dx;
                Y = voxelxyz(2) + dy;
                Z = voxelxyz(3) + dz;
                if X > imageSize(1); continue; end
                if Y > imageSize(2); continue; end
                if Z > imageSize(3); continue; end
                if X <= 0; continue; end
                if Y <= 0; continue; end
                if Z <= 0; continue; end
                cVoxelXYZ = [X; Y; Z; 1];
                cWorldXYZ = Affine * cVoxelXYZ;
                Distance = sqrt(sum((cWorldXYZ(1:3) - worldxyz').^2));
                if  Distance <= Radi
                    roiSize = roiSize + 1;
                    roiVoxelXYZ(roiSize, 1:3) = round(cVoxelXYZ(1:3));
                end
            end
        end
    end
    tempStr = string(num2str(roiVoxelXYZ));
    [tempStr, tempI] = unique(tempStr);
    roiVoxelXYZ = str2num(char(tempStr));
end

