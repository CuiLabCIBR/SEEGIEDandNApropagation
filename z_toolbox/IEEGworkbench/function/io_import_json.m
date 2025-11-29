function jsonData = io_import_json(fileName)
    % fileName  filename in JSON extension
    fid = fopen(fileName); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    jsonData = jsondecode(str); % Using the jsondecode function to parse JSON from string
end