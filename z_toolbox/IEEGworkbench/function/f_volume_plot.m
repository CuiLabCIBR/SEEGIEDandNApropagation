function [H] =  f_volume_plot(volume, Slices, dim, affine)
         [x, y, z] = meshgrid(0.5:dim(1), 0.5:dim(2), 0.5:dim(3));
         VoxCoor = [x(:), y(:), z(:)];
        RefCoor = f_vox2ref(VoxCoor, affine);
        xRef = reshape(RefCoor(:, 1), size(x));
        yRef = reshape(RefCoor(:, 2), size(y));
        zRef = reshape(RefCoor(:, 3), size(z));
        H = slice(xRef, yRef, zRef, permute(volume, [2, 1,  3]), Slices(1), Slices(2), Slices(3), 'nearest');
        shading flat
        colormap([jet(100); gray(100)])
end