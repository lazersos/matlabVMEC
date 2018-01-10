function [ output_args ] = make_coilopt_plot( varargin )
%MAKE_COILOPT_PLOT([plottype])  Plots the data from a COILOPT run
%   The MAKE_COILOPT_PLOT routine makes various plots from the COILOPT++
%   code.
%   Options:
%       'bnormal':    Plots bnormal.
%   Usage:
%       make_coilopt_plot('bnormal');
%
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.0
%   Date:       01/20/15





% Initialize some variables
plottype=0;
% Handle varargin
if nargin > 0
    for i=1:nargin
        switch varargin{i}
            case 'bnormal'
                plottype=0;
            case 'bnormal_comp'
                plottype=1;
            case 'bnormal_init'
                plottype=2;
        end
    end
end

% Load data
fid=fopen('b_norm_eq.dat','r');
temp=fscanf(fid,'%d',2);
nu=temp(1);nv=temp(2);
dataeq=fscanf(fid,'%e',[3 inf]);
fclose(fid);
ueq=reshape(dataeq(1,:),[nu nv]);
veq=reshape(dataeq(2,:),[nu nv]);
dataeq=reshape(dataeq(3,:),[nu nv]);
fid=fopen('b_norm_init.dat','r');
temp=fscanf(fid,'%d',2);
nu=temp(1);nv=temp(2);
data0=fscanf(fid,'%e',[3 inf]);
fclose(fid);
u0=reshape(data0(1,:),[nu nv]);
v0=reshape(data0(2,:),[nu nv]);
data0=reshape(data0(3,:),[nu nv]);
try
    fid=fopen('b_norm_final.dat','r');
    temp=fscanf(fid,'%d',2);
    nu=temp(1);nv=temp(2);
    dataf=fscanf(fid,'%e',[3 inf]);
    fclose(fid);
    uf=reshape(dataf(1,:),[nu nv]);
    vf=reshape(dataf(2,:),[nu nv]);
    dataf=reshape(dataf(3,:),[nu nv]);
end
try
    spline=read_coilopt_spline('fd.spline');
end





switch plottype
    case{0}
        subplot(1,3,1);
        pixplot(ueq,veq,dataeq);
        set(gca,'XTick',[0 0.5 .95],'XTickLabel',{'0' '\pi' '2\pi'});
        set(gca,'YTick',[0 0.5 .95],'YTickLabel',{'0' '\pi' '2\pi'});
        title('Equilibrium B-Normal');
        ylabel('Poloidal Angle');
        xlabel('Toroidal Angle');
        hold on;
        for i=1:spline.ncoil
            fnplt(spline.sp{i},'k');
        end
        subplot(1,3,2);
        pixplot(u0,v0,data0);
        set(gca,'XTick',[0 0.5 .95],'XTickLabel',{'0' '\pi' '2\pi'});
        set(gca,'YTick',[0 0.5 .95],'YTickLabel',{'0' '\pi' '2\pi'});
        title('Initial B-Normal');
        xlabel('Toroidal Angle');
        subplot(1,3,3);
        pixplot(uf,vf,dataf);
        set(gca,'XTick',[0 0.5 .95],'XTickLabel',{'0' '\pi' '2\pi'});
        set(gca,'YTick',[0 0.5 .95],'YTickLabel',{'0' '\pi' '2\pi'});
        title('Final B-Normal');
        xlabel('Toroidal Angle');
    case{1}
        subplot(2,2,1);
        pixplot(data0);
        subplot(2,2,2);
        pixplot(dataf);
        subplot(2,2,3);
        pixplot(dataeq-data0);
        subplot(2,2,4);
        pixplot(dataeq-dataf);
    case{2}
        subplot(1,2,1);
        pixplot(ueq,veq,dataeq);
        set(gca,'XTick',[0 0.5 .95],'XTickLabel',{'0' '\pi' '2\pi'});
        set(gca,'YTick',[0 0.5 .95],'YTickLabel',{'0' '\pi' '2\pi'});
        title('Equilibrium B-Normal');
        ylabel('Poloidal Angle');
        xlabel('Toroidal Angle');
        hold on;
        for i=1:ncoil0
            plot(spline{i}.nodes(1,:),spline{i}.nodes(2,:),'ko-');
        end
        c_save=caxis;
        subplot(1,2,2);
        pixplot(u0,v0,data0);
        set(gca,'XTick',[0 0.5 .95],'XTickLabel',{'0' '\pi' '2\pi'});
        set(gca,'YTick',[0 0.5 .95],'YTickLabel',{'0' '\pi' '2\pi'});
        title('Initial B-Normal');
        xlabel('Toroidal Angle');
        caxis(c_save);
end
        