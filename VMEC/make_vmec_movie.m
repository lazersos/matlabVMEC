function make_vmec_movie(filename,varargin)
% MAKE_VMEC_MOVIE(filename[,'2d','3d']) Create VMEC movie from file.
% 
% This function creates a movie of a VMEC equilibrium run from a VMEC movie
% file.  The plot can either be a 2D plot of the flux surfaces or a 3D plot
% of the outer most flux surface.
%
% Inputs
%   filename:  Filename of the VMEC movie file.
%
% Optional Arguments
%   'single':   2D plot of the phi=0 flux surface.
%   '2d':       2D plot of the flux surfaces phi=0,fp/4,fp/2. (default)
%   '3d':       3D plot of the outer most flux surface.
%
% Exmaple Usage
%   make_vmec_movie('movie.test','3d');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% This section handles optional input arguments.
warning off  MATLAB:divideByZero
plottype='2d';
%DIIID
xlim1=0.25; xlim2=2.75;
ylim1=-1.25; ylim2=1.25;
% LHD
%xlim1=3.0; xlim2=4.5;
%ylim1=-1.5; ylim2=1.5;
% MAST
%xlim1=0.0; xlim2=2.0;
%ylim1=-1.5; ylim2=1.5;
% LTX
%xlim1=0.0; xlim2=1.2;
%ylim1=-0.6; ylim2=0.6;
% NCSX
%xlim1=1.25; xlim2=2.75;
%ylim1=-0.75; ylim2=0.75;
% NSTX
%xlim1=-1.0; xlim2=2.5;
%ylim1=-1.75; ylim2=1.75;
% HSX
%xlim1=0.5; xlim2=1.8;
%ylim1=-0.65; ylim2=0.65;
% W7X
%xlim1=4.5; xlim2=7.0;
%ylim1=-1.5; ylim2=1.5;
%
%xlim1=0.0; xlim2 = 2.5;
%ylim1=-1.25; ylim2 = 1.25;
%ITER
%xlim1=4.0; xlim2=9;
%ylim1=-5; ylim2=5;
%BIG_TOK
%xlim1=99.0; xlim2=101.0;
%ylim1=-1.0; ylim2=1.0;
%MED_TOK
%xlim1=9.0; xlim2=11.0;
%ylim1=-1.0; ylim2=1.0;
%SMALL_TOK
%xlim1=1.0; xlim2=3.0;
%ylim1=-1.0; ylim2=1.0;
%A3_TOK
%xlim1=2.0; xlim2=4.0;
%ylim1=-1.0; ylim2=1.0;
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case 'single'
                plottype='single';
            case '2d'
                plottype='2d';
            case '3d'
                plottype='3d';
            case 'PERT_CYL'
                plottype='PERT_CYL';
            case 'MODES'
                plottype='MODES';
        end
    end
end
% Open the file
fid=fopen(filename);
% Open the movie
movname=filename(strfind(filename,'.')+1:size(filename,2));
movname=[movname '.avi'];
mov=VideoWriter(movname);
mov.Quality= 100;
open(mov);
% Size references
% 1080 1920x1080
% 720 1280x720
% 480 640x480
fig=figure('DoubleBuffer','on','Position',[1 1 1080 1080],...
        'Color','black','BackingStore','on','MenuBar','none',...
        'Name',filename,'InvertHardcopy','off');
colormap jet;
% Plot the data
theta=0:2*pi/59:2*pi;
zeta=0:2*pi/59.:2*pi;
firstframe=1;
while ~feof(fid)
    %Read file
    outtype=fscanf(fid,'%d',1);
    if isempty(outtype)
        outtype=99;
    end
    switch outtype
        case 0    %Regular data
            iter=fscanf(fid,'%d',1);
            for i=1:ns
                data=fscanf(fid,'%f',[mnmax 6]);
                rmnc(:,i)=data(:,1);
                rmns(:,i)=data(:,2);
                zmnc(:,i)=data(:,3);
                zmns(:,i)=data(:,4);
                lmnc(:,i)=data(:,5);
                lmns(:,i)=data(:,6);
            end
        case 1    %Header data
            data=fscanf(fid,'%d',[6]);
            nfp=data(1);
            ns=data(2);
            mpol=data(3);
            ntor=data(4);
            mnmax=data(5);
            niter=data(6);
            data=fscanf(fid,'%f',2*mnmax);
            xm=data(1:mnmax)';
            xn=data(mnmax+1:2*mnmax)';
            rmnc=zeros(mnmax,ns);
            rmns=zeros(mnmax,ns);
            zmnc=zeros(mnmax,ns);
            zmns=zeros(mnmax,ns);
            lmnc=zeros(mnmax,ns);
            lmns=zeros(mnmax,ns);
            continue
            %fscanf(fid,'%d',1);
            %iter=fscanf(fid,'%d',1);
            %for i=1:ns
            %    data=fscanf(fid,'%f',[mnmax 6]);
            %    rmnc(:,i)=data(:,1);
            %    rmns(:,i)=data(:,2);
            %    zmnc(:,i)=data(:,3);
            %    zmns(:,i)=data(:,4);
            %    lmnc(:,i)=data(:,5);
            %    lmns(:,i)=data(:,6);
            %end
        case 2    %New radial grid
            data=fscanf(fid,'%d',1);
            ns=data;
            data=fscanf(fid,'%f',2*mnmax);
            xm=data(1:mnmax)';
            xn=data(mnmax+1:size(data,1))';
            rmnc=zeros(mnmax,ns);
            rmns=zeros(mnmax,ns);
            zmnc=zeros(mnmax,ns);
            zmns=zeros(mnmax,ns);
            lmnc=zeros(mnmax,ns);
            lmns=zeros(mnmax,ns);
            continue
            %line1=fgetl(fid);
            %line2=fgetl(fid);
            %data=sscanf(line2,'%d');
            %iter=data;
            %for i=1:ns
            %    data=fscanf(fid,'%f',[mnmax 6]);
            %    rmnc(:,i)=data(:,1);
            %    rmns(:,i)=data(:,2);
            %    zmnc(:,i)=data(:,3);
            %    zmns(:,i)=data(:,4);
            %    lmnc(:,i)=data(:,5);
            %    lmns(:,i)=data(:,6);
            %end
        otherwise
            disp('ERROR:  Unknown line type!');
            disp(['      ' num2str(outtype)]);
            return
    end
    %Do the Fourier transformations to real space
    switch plottype
        case '2d'
            zeta=[0 pi/nfp/2 pi/nfp];
        case '3d'
            zeta=0:2*pi/(71):2*pi;
        case 'single'
            %zeta=pi/2;
            zeta=0;
        case 'PERT_CYL'
            zeta=0;
            set(fig,'Position',[1 1 640 480]);
    end
    switch plottype
        case '2d'
            theta=0:2*pi/(89):2*pi;
            r=cfunct(theta,zeta,rmnc,xm,xn);
            z=sfunct(theta,zeta,zmns,xm,xn);
            r=r+sfunct(theta,zeta,rmns,xm,xn);
            z=z+cfunct(theta,zeta,zmnc,xm,xn);
            cla;
            %Make plot 2D
            subplot(1,3,1);
            for i=1:3
                %Reset Current Axis
                subplot(1,3,i);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                plot(squeeze(r(ns,:,i)),squeeze(z(ns,:,i)),'Color','red',...
                    'LineWidth',4)
                hold on
                set(gca,'Color','black');
                axis equal
                xlim([xlim1 xlim2]);
                ylim([ylim1 ylim2]);
                parfor j=2:ns-1
                    plot(squeeze(r(j,:,i)),squeeze(z(j,:,i)),'Color','white',...
                        'LineWidth',4)
                end
                plot(squeeze(r(1,:,i)),squeeze(z(1,:,i)),'+','Color','white',...
                    'LineWidth',4,'MarkerSize',10);
                deltax=xlim1*.1;
                deltay=ylim2*.1;
                if (i==1),text(xlim1+deltax,ylim2-deltay,'0 Field Period','Color','white','FontSize',20);end
                if (i==2),text(xlim1+deltax,ylim2-deltay,'1/4 Field Period','Color','white','FontSize',20);end
                if (i==3),text(xlim1+deltax,ylim2-deltay,'1/2 Field Period','Color','white','FontSize',20);end
                xlabel('Radial Distance (R)');
                if (i==2),title(['VMEC Flux Surfaces (iteration: ' num2str(iter) ')'],'Color','white');end
                ylabel('Vertical Distance (Z)');
                hold off
            end
        case '3d'
            theta=0:2*pi/(89):2*pi;
            r=cfunct(theta,zeta,rmnc,xm,xn);
            z=sfunct(theta,zeta,zmns,xm,xn);
            r=r+sfunct(theta,zeta,rmns,xm,xn);
            z=z+cfunct(theta,zeta,zmnc,xm,xn);
            cla;
            %Make plot 3D
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            isotoro(r,z,-zeta,ns);
            xlim([-xlim2 xlim2]);
            ylim([-xlim2 ylim2]);
            zlim([ylim1 ylim2]);
            title(['VMEC Flux Surfaces (iteration: ' num2str(iter) ')'],'Color','white');
        case 'single'
            theta=0:2*pi/(89):2*pi;
            r=cfunct(theta,zeta,rmnc,xm,xn);
            z=sfunct(theta,zeta,zmns,xm,xn);
            r=r+sfunct(theta,zeta,rmns,xm,xn);
            z=z+cfunct(theta,zeta,zmnc,xm,xn);
            cla;
            % Make single 2D plot
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            plot(squeeze(r(ns,:)),squeeze(z(ns,:)),'Color','red',...
                'LineWidth',0.5)
            hold on
            set(gca,'Color','black');
            axis equal
            xlim([xlim1 xlim2]);
            ylim([ylim1 ylim2]);
            for j=2:ns-1
                plot(squeeze(r(j,:)),squeeze(z(j,:)),'Color','white',...
                    'LineWidth',0.5)
            end
            plot(squeeze(r(1,:)),squeeze(z(1,:)),'+','Color','white',...
                'LineWidth',0.5,'MarkerSize',10);
            xlabel('Radial Distance (R)');
            title(['VMEC Flux Surfaces (iteration: ' num2str(iter) ',ns=' num2str(ns) ')'],'Color','white');
            ylabel('Vertical Distance (Z)');
            hold off
        case 'PERT_CYL'
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            dr=0.5.*(rmnc(11,:)+rmnc(5,:)+zmns(11,:)-zmns(5,:));
            plot(0:1/(ns-1):1,dr./(max(dr)),'o','Color','white','LineWidth',2);
            title(['VMEC Cylindrical Perturbation (iteration: ' num2str(iter) ')'],'Color','white');
            text(0.8,0.1,['NS = ' num2str(ns)],'Color','white','FontSize',20);
            ylim([0 1.0]);
            xlim([0 1]);
            axis square;
            xlabel('Normalized Toroidal Flux')
            ylabel('Normalized Perturbation');
        case 'MODES'
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            plot(0:1/(ns-1):1,rmnc'./2e-2);
            title(['VMEC Cylindrical Perturbation (iteration: ' num2str(iter) ')'],'Color','white');
            text(0.8,0.1,['NS = ' num2str(ns)],'Color','white','FontSize',20);
            ylim([-1.5 1.5]);
            xlim([0 1]);
            axis square;
            xlabel('Normalized Toroidal Flux')
            ylabel('Normalized Perturbation');
    end
    %Add frame to movie
    if firstframe
        pause(1.0);
        firstframe=0;
    end
    %pause(.01);
    %frame=getframe(fig);
    pause(.01);
    %writeVideo(mov,frame);
end
% Close the movie
close(mov);
% Close the file
fclose(fid);
% Make a poster frame
saveas(gca,strcat(movname,'_poster.fig'),'fig');
hold on
text(0,0,'DONE','Color','yellow','FontSize',50)
hold off
clear
warning on  MATLAB:divideByZero
end
