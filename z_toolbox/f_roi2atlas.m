function [roiAtlas, p, roiAtlasIdx] = f_roi2atlas(roiVoxelXYZ, atlasImage, atlasValues, atlasLabels)
    for nv = 1:size(roiVoxelXYZ, 1)
        x = roiVoxelXYZ(nv, 1);
        y = roiVoxelXYZ(nv, 2);
        z = roiVoxelXYZ(nv, 3);
        value(nv) = atlasImage.anatomy(x, y, z);
    end
    value(value==0) = [];
    if ~isempty(value)
        valueC = tabulate(value);
        [~, i] = max(valueC(:, 3));
        j = find(atlasValues==valueC(i, 1));
        roiAtlasIdx = valueC(i, 1);
        roiAtlas = atlasLabels{j};
        p = valueC(i, 2)/size(roiVoxelXYZ, 1);
    else
        roiAtlasIdx = 0;
        roiAtlas = 'none';
        p = 1;
    end
end