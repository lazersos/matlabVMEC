function data2=optimize_coil(data,nnew)
%OPTIMIZE_COIL(data) Optimizes coils in a coil structure
%   This function takes a coil (fluxloop) data structure and attempts to
%   optimize the representation of the coil using nnew points along the
%   trajectory of the coil.  The data structure which is returned contains
%   the new coils.
nloops = data.nloops;
data2  = data;
for i=1:nloops
    npts = size(data.loops{i},2);
    x = data.loops{i}(1,:);
    y = data.loops{i}(2,:);
    z = data.loops{i}(3,:);
    dx = diff(x);
    dy = diff(y);
    dz = diff(z);
    % Assuming dl equidistant seems to result in better shapes than
    % calculing the proper dl according to length along the coil
    %ds = sqrt(dx.*dx+dy.*dy+dz.*dz);
    %dl = [0 cumsum(ds)];
    %dl = dl./max(dl);
    dl = 0:1./(npts-1):1;
    x_spl = pchip(dl,x);
    y_spl = pchip(dl,y);
    z_spl = pchip(dl,z);
    x0 = 0:1/(nnew-1):1;
    if (max(abs(diff(z))) > 0) % Avoid axisymmetric coils
        [xfinal, fval, ier] = fminsearch(@(xa) fcn_length(xa,x,y,z,x_spl,y_spl,z_spl,dl),x0(2:length(x0)-1));
        xfinal = sort(xfinal);
        xfinal = [0 xfinal 1.0];
        disp([' - Loop: ' data.loopname{i}' '  fval:' num2str(fval,'%10.2f') ' ier:' num2str(ier,'%2d')]);
    else
        xfinal = x0;
        disp([' - Skipping Loop: ' data.loopname{i}']);
    end
    hold on;
    plot3(x,y,z);
    x_new = ppval(x_spl,xfinal);
    y_new = ppval(y_spl,xfinal);
    z_new = ppval(z_spl,xfinal);
    temp(1,:) = x_new;
    temp(2,:) = y_new;
    temp(3,:) = z_new;
    data2.loops{i}=temp;
    data2.nels(i) = size(temp,2);
    plot3(x_new,y_new,z_new,'o-r');
    hold off;
    pause(0.01);
end
n_old = sum(data.nels);
n_new = sum(data2.nels);
disp(['  -  % Difference: ' num2str(100.*n_new/n_old,'%6.2f') ]);
end

function b = fcn_length(x0,x,y,z,x_spl,y_spl,z_spl,dl)
   type = 'linear';
   x0 = mod(x0,1.0);
   x0 = sort(x0);
   i = 1;
   while i< length(x0)
       if (x0(i+1) == x0(i))
           x0 = [x0(1:i) x0(i+2:length(x0))];
       end
       i = i + 1;
   end
   x_temp = ppval(x_spl,[0 x0 1]);
   y_temp = ppval(y_spl,[0 x0 1]);
   z_temp = ppval(z_spl,[0 x0 1]);
   x_new  = interp1([0 x0 1],x_temp,dl,type);
   y_new  = interp1([0 x0 1],y_temp,dl,type);
   z_new  = interp1([0 x0 1],z_temp,dl,type);
   b      = abs(sum(x-x_new) + sum(y-y_new) + sum(z-z_new));
   %b      = b + abs(sum(diff(x0)-1./length(x0+1))); % Not bad but not great either
end