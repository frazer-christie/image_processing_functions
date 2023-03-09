function nc2geotiff(varargin) 
% nc2geotiff inputs spatial information read in from NetCDF and saves data
% as a GeoTIFF (TIFF and accompaying .TFW file). This function requires the 
% use of the function imwrite2tif available at: 
% https://uk.mathworks.com/matlabcentral/fileexchange/30519-export-image-to-tif-or-tiff-file-of-selected-data-type
%
%% Syntax
% 
% nc2geotiff(imgdata,filename,x,y,datatype)
%
%
%% Description 
% 
% nc2geotiff exports NetCDF data as a GeoTIFF file and accompaying TIFF 
% world file (.tfw). It has been tested using data provided in EPSG:3031 
% Antarctic Polar Stereographic projection in NetCDF format 
% (e.g., https://nsidc.org/data/nsidc-0756/versions/3), and requires
% the following input parameters:
% 
% 'IMGDATA' is the numeric array (ususally of the form m*n) to be exported 
% 
% 'FILENAME' is a character string of the file output name. Note: no need 
% to append '.tif'.
% 
% 'X' is a numeric array of the grid dimensions of imgdata in the x 
% direction
% 
% 'Y' is a numeric array of the grid dimensions of imgdata in the y 
% direction
% 
% 'DATATYPE' is a character string specifying the data type for export. 
% Supported data types include 'logical', 'uint8', 'int8', 'uint16', 
% 'int16', 'uint32','int32', 'uint64', 'int64','single' and 'double'.
% 
%% Note:
%       Overwriting of the existing image files is not checked. Be cautious
%       with the export image file name.
% 
%% Example: 
% Read in BedMachine v3 bed data and grid information and export double 
% precision GeoTIFF with the name 'bed.tif'
%
           % imgdata = ncread('BedMachineAntarctica-v3.nc','bed')';
           % x = double(ncread('BedMachineAntarctica-v3.nc','x')); 
           % y = double(ncread('BedMachineAntarctica-v3.nc','y')); 
           % nc2geotiff(imgdata,'bed',x,y,'double')
%
%% Author Info
% 
% Code written by Frazer D.W. Christie, Scott Polar Research Institute
%
%

%% Function
% Check number of input parameters
if nargin<5
    error('Invalid number of input arguments.');
end

% assign input argument
imgdata  = varargin{1};
filename   = varargin{2};
x   = varargin{3};
y = varargin{4};
datatype = varargin{5};

% check errors
if ~isnumeric(imgdata)
     error('The first input argument (image data) must be a numeric array.');
end
if ~ischar(filename)
    error('The second input argument (filename) must be a character string.');
end
if ~isnumeric(x)
    error('The third input argument (grid x dimensions) must be a numeric array.');
end
if ~isnumeric(y)
    error('The fourth input argument (grid y dimensions) must be a numeric array.');
end
for n=1:size(imgdata,3)
    switch lower(datatype)
    case 'logical'
        BitsPerSample = 1;
        SampleFormat  = 1;
        imgdata = logical(imgdata);
    case 'uint8'
        BitsPerSample = 8;   
        SampleFormat  = 1;         
        imgdata = uint8(imgdata);    
    case 'int8'
        BitsPerSample = 8;
        SampleFormat  = 2;       
        imgdata = int8(imgdata);          
    case 'uint16'
        BitsPerSample = 16;       
        SampleFormat  = 1;         
        imgdata = uint16(imgdata);           
    case 'int16'
        BitsPerSample = 16;       
        SampleFormat  = 2;         
        imgdata = int16(imgdata);        
    case 'uint32'
        BitsPerSample = 32;       
        SampleFormat  = 1;                
        imgdata = uint32(imgdata);            
    case 'int32'
        BitsPerSample = 32;       
        SampleFormat  = 2;     
        imgdata = int32(imgdata);         
     case 'single'
        BitsPerSample = 32;        
        SampleFormat  = 3;                
        imgdata = single(imgdata);       
    case {'uint64','int64','double'}
        BitsPerSample = 64;        
        SampleFormat  = 3;                
        imgdata = double(imgdata);           
    otherwise
        error('Invalid output data type.');
    end

% Write GeoTIFF image
imwrite2tif(imgdata,[],strcat(filename,'.tif'),datatype); % Note: imwrite2tif is used becuase MATLABs default geotiffwrite is rather restrictive, and truncates all double precision data to uint8 (1-255).

% Assemble accompanying .tfw file using grid x and y information provided 
% in the NetCDF file. 
xlim = [min(x(:)),max(x(:))]; % Calculate min/max x grid limits 
ylim = [min(y(:)),max(y(:))]; % Calculate min/max y grid limits 

Atfw = round((max(x)-min(x))./(size(imgdata,2)-1)); % cellsize in x-direction  
Dtfw = 0; %rotation about y-axis
Btfw = 0; %rotation about x-axis
Etfw = round((max(y)-min(y))./(size(imgdata,1)-1).*-1); % cellsize in y-direction (almost always negative)
Ctfw = min(x); % center x-coordinate of upper left pixel
Ftfw = max(y); % center y-coordinate of upper left pixel

tfw_data = vertcat(Atfw,Dtfw,Btfw,Etfw,Ctfw,Ftfw); % concatenate TFW contents
writematrix(tfw_data,strcat(filename,'.txt')); % begin by making .txt file

directory = pwd;
filelist = dir(strcat(filename,'.txt')) % find .txt file in directory

% change extension from .txt to .tfw
for i = 1:numel(filelist)
    file = fullfile(directory,filelist(i).name);
    [tempDir,tempFile] = fileparts(file);
    status = copyfile(file, fullfile(tempDir,[tempFile,'.tfw']));  
end 

txt = strcat(filename,'.txt'); % Get hamdle on now uneeded .txt file
delete(txt) % Delete now uneeded .txt file

end
