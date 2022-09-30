function K = getEVs_loren(nEV, ind, solarFileHdr, solarFile)
    % Import_script
    
    % Replace previous call to Import_script.m with flexible code
    hdr_data = read_envihdr(solarFileHdr);

    new_data=multibandread(solarFile,[hdr_data.lines,hdr_data.samples,hdr_data.bands],...
        strcat(hdr_data.format,'=>float32'),hdr_data.header_offset,...
        hdr_data.interleave,hdr_data.machine);

    % get rid of 3rd dimension (resulting 'spectra' object should be n x 2160 double):
    if size(new_data, 1) == 500 && size(new_data,2)== 1
        spectra(:,:) = double(new_data(:,1,:));
    elseif size(new_data, 1) == 2160
        spectra = new_data';
    end

    FR_index = ind;
    [U,s,V] = svd(spectra(:,FR_index),'econ');
    EV =  V(:,1:nEV);
    xx = (ind-mean(ind))/100;
    % Just using constant SIF shape here...
    K = [EV V(:,1).*xx.^1' V(:,1).*xx.^2' V(:,1).*xx.^3' ones(length(ind),1)]; 
    % Note: in Matlab, A.*B is the element-by-element product of A and B
    % and A.^B is the matrix with elements A(i,j) to the B(i,j) power.
end