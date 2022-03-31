function data = read_travis
%read_travis Reads a TRAVIS directory
%   The READ_TRAVIS routine reads the files in the current directory into a
%   a structure for plotting.

files=dir('beamtrace*');
for i=1:length(files)
    temp = importdata(files(i).name,' ',1);
    temp = temp.data;
    for j=1:max(temp(:,1)) % Loop over rays
        dex = temp(:,1)==j;
        data.path{j,i} = temp(dex,2);
        data.x{j,i} = temp(dex,3);
        data.y{j,i} = temp(dex,4);
        data.z{j,i} = temp(dex,5);
        data.nx{j,i} = temp(dex,6);
        data.ny{j,i} = temp(dex,7);
        data.nz{j,i} = temp(dex,8);
        data.rho{j,i} = temp(dex,9);
        data.ne{j,i} = temp(dex,10);
        data.te{j,i} = temp(dex,11);
        data.B{j,i} = temp(dex,12);
    end
end
files=dir('Pabs_Icd_profiles*');
for i=1:length(files)
    temp = importdata(files(i).name,' ',1);
    temp = temp.data;
    data.reff(:,i) = temp(:,1);
    data.dPpdV(:,i) = temp(:,2);
    data.dPtdV(:,i) = temp(:,3);
    data.Pp(:,i) = temp(:,4);
    data.Pt(:,i) = temp(:,5);
    data.dP0dV(:,i) = temp(:,6);
    data.P0(:,i) = temp(:,7);
    data.jcd_p(:,i) = temp(:,8);
    data.jcd_t(:,i) = temp(:,9);
    data.Itor_p(:,i) = temp(:,10);
    data.Itor_t(:,i) = temp(:,11);
    data.effCD(:,i) = temp(:,12);
end

return;
end

