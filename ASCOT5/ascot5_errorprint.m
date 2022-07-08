function ascot5_errorprint(a5file,runid)
%ASCOT5_ERRORPRINT Prints the endcondition and errors from an ASCOT 5 run.
%   The ASCOT5_ERRORPRINT subroutine prints the end and  error conditions 
%   from an ASCOT5 run.
%
%   Example:
%       a5file='ascot5_W7X_20180821_012_5100_h8.h5';
%       id=0614695199;
%       ascot5_errorprint(a5file,id);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

    % Use active run
    if isempty(runid)
        try
            runid=h5readatt(a5file,'/results','active');
            disp(['  Using runid: ' runid]);
        catch
            runid=[];
            return;
        end
    end

% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

if isempty(runid)
    runid=h5readatt(a5file,'/results','active');
    disp(['  Using runid: ' runid]);
end


% from error.h
file_errs={'mccc_wiener.c','mccc_push.c','mccc_coefs.c','mccc.c','step_fo_vpa.c',...
    'step_gc_cashkarp.c','step_gc_rk4.c','N0_3D.c','N0_ST.c','B_3DS.c',...
    'B_2DS.c','B_STS.c','B_GS.c','plasma_1D.c','plasma_1DS.c',...
    'plama.c','E_field.c','neutral.c','E_1DS.c','B_field.c',...
    'particle.c'};

endpath = ['/results/run_' num2str(runid,'%10.10i') '/endstate/'];
try
    endcond = h5read(a5file,[endpath '/endcond']);
catch
    disp(['ERROR: Could not find run number or endstate: ' num2str(runid,'%10.10i')]);
    return;
end
errorline = h5read(a5file,[endpath '/errorline']);
errormod = h5read(a5file,[endpath '/errormod']);
errormsg = h5read(a5file,[endpath '/errormsg']);

end_str='';

disp('END CONDITIONS');
end_range=unique(endcond);
for i = end_range'
    if i==0, continue;end
    num_part = sum(endcond == i);
    
    % from endcond.h ASCOT5
    if (i==1)
        end_str='Maximum Simulation Time';
    elseif (i==2)
        end_str='Minimum Energy';
    elseif (i==4)
        end_str='Thermalized';
    elseif (i==8)
        end_str='Wall collision';
    elseif (i==16)
        end_str='Minimum rho';
    elseif (i==32)
        end_str='Maximum rho';
    elseif (i==64)
        end_str='Poloidal limit';
    elseif (i==128)
        end_str='Toroidal limit';
    elseif (i==256)
        end_str='Wall time exceeded';
    elseif (i==512)
        end_str='Hybrid mode condition';
    else
        end_str=['Unknown: ' num2str(i,'%i')];
    end
    
    disp(['     Detected ' num2str(num_part,'%i') ' particles with end condition: ' end_str]);
end

% Now handle errors
if any(errormsg~=0)
    num_part = sum(errormsg~=0);
    disp(['     Detected ' num2str(num_part,'%i') ' particles with errors!']);
    disp('ERRORS');
    error_range=unique(errormsg)';
    
    % from error.c
    for i = error_range
        if i==0, continue;end
        err_dex = errormsg == i;
        num_part = sum(err_dex);
        if (i==1)
            err_str='Input evaluation failed (marker could be outside input data grid).';
        elseif (i==2)
            err_str='Input was not recognized (offload or target data could uninitalized).';
        elseif (i==3)
            err_str='Input evaluation yeilds unphysical results (something could be wrong with the input).';
        elseif (i==4)
            err_str='One or more markers fields are unphysical or inconsistent (marker input could be corrupted).';
        elseif (i==5)
            err_str='Time step is zero, NaN, or smaller than MIN_ALLOWED_TIME_STEP';
        elseif (i==6)
            err_str='Wiener array is full of rejected steps or limits could be too conservative or initial step too large.';
        else
            err_str=[' Unknown error value: ' num2str(i,'%i')];
        end
        disp(['     Detected ' num2str(num_part,'%i') ' particles with error condition: ' err_str]);
        mod_range = unique(errormod(err_dex))';
        for i = mod_range
            mod_dex = and(err_dex,errormod == i);
            num_part = sum(mod_dex);
            line_range = unique(errorline(mod_dex))';
            for j = line_range
                num_part = sum(and(mod_dex,errorline == j));
                disp(['          ' num2str(num_part,'%i') ' from line ' num2str(j,'%i') ...
                    ' in ' file_errs{i}]);
            end
            
        end
    end
    
end

end