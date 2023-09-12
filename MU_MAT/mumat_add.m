function data = mumat_add(datain1,datain2)
%MUMAT_ADD Combines two mu material structures into one
%   MUMAT_ADD takes two mu material structures as read by READ_MUMAT and 
%   returns a third mu material structure which contains the elements
%   of the first two.
%
%   Example:
%       mu_data1=read_mumat('solid_object.dat');
%       mu_data2=read_mumat('solid_object_more.dat');
%       mu_data_all=mumat_add(mu_data1,mu_data2);
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           9/12/2023

data = [];
if isempty(datain1) && isempty(datain2)
    return;
elseif isempty(datain1)
    data=datain2;
    return;
elseif isempty(datain2)
    data=datain1;
    return;
end

data.machine = datain1.machine;
data.date = datestr(today);
data.coords = [datain1.coords datain2.coords];
data.tet = [datain1.tet datain2.tet+datain1.nvertex];
data.func_dex = [datain1.func_dex datain2.func_dex+datain1.nstate];
data.state_func = [datain1.state_func datain2.state_func];
return
end