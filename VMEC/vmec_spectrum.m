function tmp=vmec_spectrum(data,varargin)
%VMEC_SPECTRUM Dumps the VMEC spectrum of a wout file to the screen.
%   The VMEC_SPECTRUM function takes a VMEC wout structure as read by the
%   READ_VMEC function.  It returns to screen the boundary harmonic
%   information in multiple formats.  The default is the VMEC INDATA
%   format.  Options include:
%       flip:       Flip the torodial angle
%       spec:       Output in SPEC format
%       nescoil:    Output in NESCOIL format
%       focus:      Output in FOCUS format (Bn=0)
%       ns:         Specify non-boundary surface for output.
%       
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.00

lspec = 0;
lvmec = 1;
lnescoil = 0;
lfocus = 0;
lhardcode = 0; % For hardcoding things in matlab
lflip = 0; % Flip the jacobian sign
s=[];

% Handle varargin
if nargin > 1
    i=1;
    while i < nargin
        switch varargin{i}
            case{'flip','FLIP'}
                lflip = 1;
            case{'spec','SPEC'}
                lspec = 1;
            case{'necoil','NESCOIL'}
                lnescoil = 1;
            case{'focus','FOCUS'}
                lfocus = 1;
            case{'ns','NS'}
                i=i+1;
                s=varargin{i};
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end
if (lspec || lnescoil || lfocus), lvmec=0; end
if (isempty(s)), s = data.ns; end

% Flip toroidal direction
if (lflip)
    for i=1:data.mnmax
        m=data.xm(i);
        if (m>0)
            if mod(m,2) == 0
                data.zmns(i,:) = -data.zmns(i,:);
            else
                data.rmnc(i,:) = -data.rmnc(i,:);
            end
        end
    end
    data.xn = -data.xn;
end

% Note that the VMEC Kernel in MATLAB is (mu+nv) based on read_vmec
if (data.iasym == 0)
    if lspec
        disp('mn   m   n   rmnc   zmns');
        for i=1:data.mnmax
            disp([num2str(i,'%3.3d') '   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.xn(i),'%+2.2d') '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') ]);
        end
        disp('  m   n   rmnc   zmns   rmns   zmnc');
        for i=1:data.mnmax
            disp(['   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.xn(i)./data.nfp,'%+2.2d') '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') '      0.0000000000E+00     0.0000000000E+00' ]);
        end
    elseif lnescoil
        disp('  m   n   rmnc   zmns   lmns   rmns   zmnc   lmnc');
        for i=1:data.mnmax
            disp(['   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.xn(i)./data.nfp,'%+2.2d') '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') '   ' num2str(data.lmns(i,s),'%e') '      0.0000000000E+00     0.0000000000E+00     0.0000000000E+00' ]);
        end
    elseif lvmec 
        rax = data.rmnc(1:data.ntor+1,1);
        zax = -data.zmns(1:data.ntor+1,1);
        disp(['  RAXIS = ' num2str(rax',' %20.12E ')]);
        disp(['  ZAXIS = ' num2str(zax',' %20.12E ')]);
        for i=1:data.mnmax
            disp(['  RBC(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.rmnc(i,s),'%20.12E') '  ZBS(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.zmns(i,s),'%20.12E') ]);
            tmp(i,:) = [-data.xn(i)./data.nfp data.xm(i) data.rmnc(i,s) data.zmns(i,s)];
        end

    elseif lfocus
        disp('#bmn  bNfp nbf');
        disp([num2str(data.mnmax,'%3i ') ' ' num2str(data.nfp,'%3i ') ' ' num2str(data.mnmax,'%3i ')]);
        disp('#------plasma boundary harmonics-------');
        disp('# n m Rbc Rbs Zbc Zbs');
        for i=1:data.mnmax
            disp(['   '  num2str(data.xn(i)./data.nfp,'%+2.2d') '   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.rmnc(i,s),'%e') '      0.0000000000E+00      0.0000000000E+00 ' num2str(-data.zmns(i,s),'%e') ]);
        end
    elseif lhardcode
        disp(['  RAXIS = [' num2str(data.rmnc(1:data.ntor+1,1)',' %20.12E ') '];']);
        disp(['  ZAXIS = [' num2str(data.zmns(1:data.ntor+1,1)',' %20.12E ') '];']);
        nmax = max(data.xn)./data.nfp+1;
        for i=1:data.mnmax
            disp(['  RBC(' num2str(nmax-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i)+1,'%2.2d') ') = ' num2str(data.rmnc(i,s),'%20.12E') ';  ZBS(' num2str(nmax-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i)+1,'%2.2d') ') = ' num2str(data.zmns(i,s),'%20.12E') ';']);
        end
    end
else
    if lspec
        disp('mn   m   n   rmnc   rmns   zmnc   zmns');
        for i=1:data.mnmax
            disp([num2str(i,'%3.3d') '   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.xn(i),'%+2.2d')...]
                '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(data.rmns(i,s),'%e')...
                '   ' num2str(data.zmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') ]);
        end
    elseif lnescoil
        disp('  m   n   rmnc   zmns   lmns   rmns   zmnc   lmnc');
        for i=1:data.mnmax
            disp(['   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.xn(i),'%+2.2d')...]
                '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') '   ' num2str(data.lmns(i,s),'%e')...
                '   ' num2str(data.zmnc(i,s),'%e') '   ' num2str(data.zmns(i,s),'%e') '   ' num2str(data.lmnc(i,s),'%e') ]);
        end
    elseif lfocus
        disp('#Nbmn  Nfp Nbnorm:');
        disp([num2str(data.mnmax,'%3i ') num2str(data.nfp,'%3i ') num2str(data.mnmax,'%3i ')]);
        disp('#------plasma boundary harmonics-------');
        disp('# n m Rbc Rbs Zbc Zbs');
        for i=1:data.mnmax
            disp(['   '  num2str(data.xn(i),'%+2.2d') '   ' num2str(data.xm(i),'%2.2d') '   ' num2str(data.rmnc(i,s),'%e') '   ' num2str(-data.rmns(i,s),'%e') '   ' num2str(data.zmnc(i,s),'%e') '   ' num2str(-data.zmns(i,s),'%e') ]);
        end
    elseif lvmec
        rax_cc = data.rmnc(1:data.ntor+1,1);
        zax_cs = -data.zmns(1:data.ntor+1,1);
        rax_cs = -data.rmns(1:data.ntor+1,1);
        zax_cc = data.zmnc(1:data.ntor+1,1);
        disp(['  RAXIS = ' num2str(rax',' %20.12E ')]);
        disp(['  ZAXIS = ' num2str(zax',' %20.12E ')]);
        disp(['  RAXIS_CC = ' num2str(rax_cc',' %20.12E %20.12E %20.12E %20.12E %20.12E ')]);
        disp(['  RAXIS_CS = ' num2str(rax_cs',' %20.12E %20.12E %20.12E %20.12E %20.12E ')]);
        disp(['  ZAXIS_CS = ' num2str(zax_cs',' %20.12E %20.12E %20.12E %20.12E %20.12E ')]);
        disp(['  ZAXIS_CC = ' num2str(zax_cc',' %20.12E %20.12E %20.12E %20.12E %20.12E ')]);
        for i=1:data.mnmax
            disp(['  RBC(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.rmnc(i,s),'%20.12E') '  ZBS(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.zmns(i,s),'%20.12E') ]);
            disp(['    RBS(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.rmns(i,s),'%20.12E') '    ZBC(' num2str(-data.xn(i)./data.nfp,'%2.2d') ',' num2str(data.xm(i),'%2.2d') ') = ' num2str(data.zmnc(i,s),'%20.12E') ]);
        end
    end
end

end

