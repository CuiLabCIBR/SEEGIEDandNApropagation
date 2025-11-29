function f_densityScatter(X, Y, varargin)
    %% Set default values
    default_xbin = 50;
    default_ybin = default_xbin;
    light_darkGreen=[0.9400    0.9700    0.9600 % Lightest green
                    0.8900    0.9300    0.9200
                    0.8200    0.9100    0.8800
                    0.6900    0.8500    0.7700
                    0.5900    0.7800    0.6900
                    0.5500    0.7500    0.6500
                    0.4500    0.6500    0.5600
                    0.4000    0.5800    0.4900
                    0.3500    0.5100    0.4200
                    0.2500    0.3600    0.3100
                    0.1300    0.1700    0.1400]; % Darkest green
    darkBlue_yellowGreen=[0.2700         0    0.3300
                        0.2700    0.2300    0.5100
                        0.1900    0.4100    0.5600
                        0.1200    0.5600    0.5500
                        0.2100    0.7200    0.4700
                        0.5600    0.8400    0.2700
                        0.9900    0.9100    0.1300];
    light_darkRed=[1.0  0.95 0.95
                      1.0  0.7  0.7
                      0.9  0.4  0.4
                      0.6  0.1  0.1
                      0.2  0.0  0.0];
    gold_purple = [0.98  0.94  0.70   % 浅金色
                      0.85  0.65  0.85   % 紫罗兰
                      0.55  0.20  0.70   % 深紫色
                      0.25  0.05  0.35]; % 近乎黑色的深紫
    black_red_yellow = [0.00  0.00  0.00   % 纯黑色
              0.20  0.00  0.00   % 暗红色
              0.50  0.00  0.00   % 深红色
              0.80  0.00  0.00   % 纯红色
              1.00  0.20  0.00   % 红橙色
              1.00  0.50  0.00   % 橙色
              1.00  0.75  0.20   % 金橙色
              1.00  0.90  0.50   % 浅金黄
              1.00  0.95  0.70   % 亮金黄
              1.00  1.00  0.90]; % 淡金黄（接近白色）
    ocean = [0.94  0.95  0.97   % 浪花白（Foam White）
              0.85  0.92  0.97   % 浅天蓝（Light Sky Blue）
              0.70  0.85  0.95   % 天蓝色（Sky Blue）
              0.55  0.75  0.90   % 浅海蓝（Light Ocean Blue）
              0.40  0.65  0.85   % 海洋蓝（Ocean Blue）
              0.25  0.50  0.75   % 中深海蓝（Medium Deep Blue）
              0.15  0.35  0.60   % 深海蓝（Deep Blue）
              0.08  0.20  0.40   % 暗深海蓝（Dark Deep Blue）
              0.02  0.08  0.25]; % 深海墨蓝（Abyssal Blue）
    colorList = {light_darkGreen, darkBlue_yellowGreen, ...
                light_darkRed, gold_purple, ...
                black_red_yellow, ocean};
    %% Parse optional inputs
    switch nargin
        case 2
            xbin = default_xbin;
            ybin = default_ybin;
            colorList = colorList{1};
        case 3
            colorList = colorList{varargin{1}};
            xbin = default_xbin;
            ybin = default_ybin;
        case 4
            colorList = colorList{varargin{1}};
            xbin = varargin{2};
            ybin = xbin;
        case 5  
            colorList = colorList{varargin{1}};
            xbin = varargin{2};
            ybin = varargin{3};
        case 6  
            xbin = varargin{2};
            ybin = varargin{3};
            colorList = varargin{4};
        otherwise
            error('Invalid number of inputs. Requires at least 2 inputs (X and Y)');
    end

    %% Compute 2D kernel density estimate
    [densityMesh, Xmesh, Ymesh] = f_ksdensity2d(X, Y, xbin, ybin);

    % Interpolate density at each original (X,Y) point
    density = interp2(Xmesh, Ymesh, densityMesh, X, Y);

    % Create color mapping function
    colorFunc=f_colorFuncFactory(colorList);

    % Normalize densities to [0,1] range and map to colors
    CData = colorFunc((density-min(density))./(max(density)-min(density)));
    
    % Create density-colored scatter plot
    scatter(X, Y, 'filled', 'CData', CData);
end