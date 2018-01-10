function vmec_data = flip_spectrum( vmec_data )
%FLIP_SPECTRUM Flips the deffinition of the poloidal angle
%   Detailed explanation goes here

vmec_data.zmns = -vmec_data.zmns;
vmec_data.bsubsmns = -vmec_data.bsubsmns;
vmec_data.xn   = -vmec_data.xn;
if (vmec_data.iasym)
    vmec_data.rmns = -vmec_data.rmns;
    vmec_data.gmns = -vmec_data.gmns;
    vmec_data.bsubumns = -vmec_data.bsubumns;
    vmec_data.bsubvmns = -vmec_data.bsubvmns;
    vmec_data.bsupumns = -vmec_data.bsupumns;
    vmec_data.bsupvmns = -vmec_data.bsupvmns;
    vmec_data.currumns = -vmec_data.currumns;
    vmec_data.currvmns = -vmec_data.currvmns;
end

end

