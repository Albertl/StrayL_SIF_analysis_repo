function SIF = fitAllFR(vegFile,solarFile)

% Make strings with .hdr extension for solar file
solarFileHdr = strcat(solarFile,'.hdr');

% Make strings with .hdr extension for veg file and load, or convert from
% table
if istable(vegFile)
    data1 = table2array(vegFile);
else
    vegFileHdr = strcat(vegFile,'.hdr');
    % load data:
    data1 = read_hdr(vegFileHdr,vegFile);
end

% Specify input file:
%file = "Spectra_for_Frankenberg/Solar_no_SL_frames.hdr"
%file = "Spectra_for_Frankenberg/Vegetation_no_SL_SIFL2.hdr"
%file = "Spectra_for_Frankenberg/Vegetation_SL_Dmat_scenario3_SIFL2.hdr"
%file = "/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Spectra_output/Vegetation_no_SL_SIFL2.hdr";
% Specify fit range:
ind = 1530:1700;
% construct K matrix:
K = getEVs_loren(3,ind,solarFileHdr,solarFile);

s = size(data1);
SIF = zeros(s(2),1);

for i=1:s(2)
    y = double(data1(ind,i));
    x = K\y;
    % Note, in Matlab, x = A\B solves the system of linear equations A*x =
    % B. see help for mldivide.
    SIF(i) = x(end);
end
end
