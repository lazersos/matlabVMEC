function poincm(rbc,zbs,vrhobc,vthebc,vphibc,cr,cthe,cphi,rmaj,nfp)
%POINCM(vr,vz,cr,cz) Poincaire plot
%   Creates a poinciare plot given a vector field

maxiter=10000;
stepsize=.01;
hold on;
%  We track each particle
for i=1:max(size(cx))
    % Start each particle at it's beginning point
    pr=cr(i);
    bigr=pr*cos(cthe(i));
    bigz=pr*sin(cthe(i));
    pthe=cthe(i);
    pzeta=cphi(i);
    partstop=0;
    n=1;
    plot(bigr,bigz,'.');
    title(['Iteration n=',num2str(n),' x(r,the,zeta)=(',num2str(pr),',',num2str(pthe),',',num2str(pzeta),')'])
    axis([4 7 -1 1]);
    pause(.1);
    while ~partstop
        % Calculate the fields at a point
        rt=cfunct(pthe,pzeta,rbc,nfp)+rmaj;
        zt=sfunct(pthe,pzeta,zbs,nfp);
        radius=sqrt(rt.*rt+zt.*zt);
        vrho=cfunct(pthe,pzeta,vrhobc,nfp);
        vthe=cfunct(pthe,pzeta,vthebc,nfp);
        vphi=cfunct(pthe,pzeta,vphibc,nfp);
        % Deal with Contravariant form
        vrho=vrho.*vphi;
        vthe=vthe.*vphi;
        % Spline to r position
        vrhoa=spline(radius,vrho,pr);
        vthea=spline(radius,vthe,pr);
        vphia=spline(radius,vphi,pr);
        % Normalize Vectors
        vnorm=sqrt(vrhoa*vrhoa+vthea*vthea+vphia*vphia);
        vrhoa=vrhoa/vnorm;
        vthea=vthea/vnorm;
        % Now update position
        dx=(vrhoa*cos(pthe)-vthea*sin(pthe));
        dy=(vthea*cos(pthe)+vrhoa*sin(pthe));
        dphi=vphia/vnorm;
        bigr=bigr+dx*stepsize;
        bigz=bigz+dy*stepsize;
        pr=sqrt(bigr*bigr+bigz*bigz);
        pthe=atan(bigz/(bigr-rmaj));
        pzeta=pzeta+dphi*stepsize;
        % Check to see if we've crossed the phi=2*pi plane
        if pzeta >= 2*pi/nfp
            plot(bigr,bigz,'.k','MarkerSize',10);
            pzeta=pzeta-2*pi/nfp;
            pause(.001);
        end
        % update number n
        n=n+1;
        % Check for exit
        if (n>maxiter)
            partstop=1;
            plot(bigr,bigz,'+k');
        end
        %if ~mod(n,10)
           title(['Iteration n=',num2str(n),' x(r,the,zeta)=(',num2str(pr),',',num2str(pthe),',',num2str(pzeta),')'])
           pause(.01);
        %end
    end
end

end

