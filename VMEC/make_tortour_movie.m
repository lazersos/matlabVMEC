function make_tortour_movie(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
author='lazerson@pppl.gov';
warning off  MATLAB:divideByZero
%read_coils(coilname,'simple');% Setup Movie File
mov=VideoWriter(filename);
mov.Quality= 100;
mov.FrameRate=30;
open(mov);
% Resolution Note
% 480 640x480 (Standard Def)
% 720 1280x720
% 1080 1920x1080
set(gcf,'DoubleBuffer','on','Position',[1 1 1920 1080],...
    'Color','black','BackingStore','on','MenuBar','none',...
    'Name',filename);
set(gca,'nextplot','replacechildren','Color','black');
axis vis3d
%  Get first frame
frame=getframe(gcf);
writeVideo(mov,frame);
pause(.1);
bpos=campos;
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[20. 0 0];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[20. 20. 0];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[0. 20. -5.];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[-20. 20. 5.];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[-20. 0. 0.];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[-20. -20. 0.];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=[0 0 30.];
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%  Now track from Start Position to X=R;
frames=90;
spos=campos;
starget=camtarget;
epos=bpos;
dx=(epos(1)-spos(1))/frames;
dy=(epos(2)-spos(2))/frames;
dz=(epos(3)-spos(3))/frames;
for i=1:frames
    camdolly(dx,dy,dz,'fixtarget','data');
    pause(.1);
    frame=getframe(gcf);
    writeVideo(mov,frame);
end
%Close the movie
close(mov);
end

