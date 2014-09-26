function [M]=PlotSphere(LR,RC,DV,VAR,Size,Color,Colorcat,Colorcatname,FIXED_max,x)
mov=0;
r=45;
PaperVersion2D=1; %Set to 0 to see the 3D version

%LR: Left-Right position
%RC: Rostrocaudal position
%DV: Dorsoventral position
%FIXED_max = 1;

%% change base coordinates
BaseShift = [1000 0 -1500];%This is the camera target if nothing is done...
%to zoom in more easily I'm changing the base so the camera target is 0 0 0
LR = LR-BaseShift(1);
RC=RC-BaseShift(2);
DV=DV-BaseShift(3);

%% Getting VAR ready for plotting.
% If VAR is always a non-null
% integer then you can use it directly in cubehelix by specifying x=1
% otherwise, the function first multiplies VARY by 10^x to make sure its
% variability is maximized after ceiling it to an integer
if nargin<10
    x=3;% here I'm doing *1000 as an example
end

%% check input variables
if nargin<8
    Colorcatname={};
end

if nargin<7
    Colorcat={};
end

if nargin<5
    Size=1;
end

if nargin<6
    Color=1;
end

%% Make sure you don't have negative values as input for cubehelix or set a
% slide value
Slide = 0;
if sum(VAR<0)>0
    Slide = -min(VAR)+ 10^(-x);
end


VAR_cube = ceil((VAR+Slide)*10^x);
%GRAD=cubehelix(max(VARZ_cube));%Note that you can choose the values of...
...your gradient by doing cubehelix_view and then specifying the values in...
    ...the function, for instance: GRAD=cubehelix(max(VARZ_cube),0.5, -0.8, 1.5, 1.7, [1,0]);
%GRAD=cubehelix(max(VAR_cube),2.56, 0, 1, 1, [1 0]);
if FIXED_max>0
    GRAD=cubehelix(FIXED_max*10^x,1.83, 0.77, 3, 1.1, [1 0]);
else
    GRAD=cubehelix(max(VAR_cube),1.83, 0.77, 3, 1.1, [1 0]);
end

%% uper bound the diameters of the spheres to 45
if Color && ~Size
    VAR_diam = VAR_cube./(max(VAR_cube)./25);
elseif Color && Size
        VAR_diam = VAR_cube./(max(VAR_cube)./45);
end

%% Plot
figure()
NU = length(LR); % this is the number of points you have in your plot
if Color && Size
    if isempty(Colorcat)
        for ss=1:NU
            [x,y,z]=sphere(r);
            sph=surf(x*VAR_diam(ss)+LR(ss), y*VAR_diam(ss)+RC(ss), z*VAR_diam(ss)+DV(ss),'EdgeColor','none','LineStyle','none','FaceColor',GRAD(VAR_cube(ss),:), 'FaceLighting', 'phong','AmbientStrength', 0.5);
            hold on
        end
    else
        %GRAD2=cubehelix(max(Colorcat),3, -2.2, 3, 0.76, [0 1]);
        if length(Colorcatname)==9
            GRAD2=[0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
        elseif length(Colorcatname)==7
            GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.8 0.8 0.8];
        elseif length(Colorcatname)==10
            %GRAD2=[0.8 0.8 0.8; 0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
            GRAD2=[0.8 0.8 0.8; 0 1 1; 0 0.8 1; 1 0.3 0; 1 0 0; 1 0 0 ; 1 0 0.6; 1 0.2 0.5; 0 1 0; 0 0 0];
        elseif length(Colorcatname)==2
            GRAD2=[0.8 0.8 0.8; 0.7 0 0.7];
        elseif length(Colorcatname)==4
            GRAD2=[0 0.5 1;0 1 1;1 0.75 0;0.5 0 0.1];
        else
            fprintf(1, 'Problem with the GRAD2 definition, the length of the color code is %d and the colormap (GRAD2) can only be 9 or 10', length(Colorcatname)); 
        end
        for ss=1:NU
            [x,y,z]=sphere(r);
            sph=surf(x*VAR_diam(ss)+LR(ss), y*VAR_diam(ss)+RC(ss), z*VAR_diam(ss)+DV(ss),'EdgeColor','none','LineStyle','none','FaceColor',GRAD2(Colorcat(ss),:), 'FaceLighting', 'phong','AmbientStrength', 0.5);
            hold on
        end
        for ii=1:length(Colorcatname)
            text(-1000,1400,1600-100.*ii,Colorcatname(ii), 'Color', GRAD2(ii,:))
        end
    end
end
if Color && ~Size
    if isempty(Colorcat)
        for ss=1:NU
            [x,y,z]=sphere(r);
            sph=surf(x*max(VAR_diam)+LR(ss), y*max(VAR_diam)+RC(ss), z*max(VAR_diam)+DV(ss),'EdgeColor','none','LineStyle','none','FaceColor',GRAD(VAR_cube(ss),:), 'FaceLighting', 'phong','AmbientStrength', 0.5);
            hold on
        end
    else
        if length(Colorcatname)==9
            GRAD2=[0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
        elseif length(Colorcatname)==7
            %GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.8 0.8 0.8];%Great with dark font
            GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.2 0.2 0.2];%Great with white font
        elseif length(Colorcatname)==10
            GRAD2=[0.8 0.8 0.8; 0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
        elseif length(Colorcatname)==2
            GRAD2=[0.8 0.8 0.8; 0.7 0 0.7];
        elseif length(Colorcatname)==4
            GRAD2=[0 0.5 1;0 1 1;1 0.75 0;0.5 0 0.1];
        else
            fprintf(1, 'Problem with the GRAD2 definition, the length of the color code is %d and the colormap (GRAD2) can only be 9 or 10', length(Colorcatname)); 
        end
        for ss=1:NU
            [x,y,z]=sphere(r);
            sph=surf(x*max(VAR_diam)+LR(ss), y*max(VAR_diam)+RC(ss), z*max(VAR_diam)+DV(ss),'EdgeColor','none','LineStyle','none','FaceColor',GRAD2(Colorcat(ss),:), 'FaceLighting', 'phong','AmbientStrength', 0.5);
            hold on
        end
        for ii=1:length(Colorcatname)
            text(-1000,1400,1600-100.*ii,Colorcatname(ii), 'Color', GRAD2(ii,:))
        end
    end
       
end

if ~Color && Size
    for ss=1:NU
        [x,y,z]=sphere(r);
        sph=surf(x*VAR_diam(ss)+LR(ss), y*VAR_diam(ss)+RC(ss), z*VAR_diam(ss)+DV(ss),'EdgeColor','none','LineStyle','none','FaceColor','y', 'FaceLighting', 'phong','AmbientStrength', 0.5);
        hold on
    end
end
    
light('Position', [-1 0 0.5], 'Style', 'infinite')

%% set the limits of the axis to see the spheres well in 3D version
ZLIM=get(gca, 'ZLim');
Zrange=ZLIM(2) - ZLIM(1);
XLIM=get(gca, 'XLim');
Xrange=XLIM(2) - XLIM(1);
YLIM=get(gca, 'YLim');
Yrange=YLIM(2) - YLIM(1);
if PaperVersion2D==1
    Maxrange=max([Zrange, Yrange]);%we don't care about the right left axis since we look at a sagittal view
else
    Maxrange=max([Zrange, Xrange, Yrange]);
    set(gca,'XLim', [XLIM(1)-(Maxrange-Xrange)/2 XLIM(2)+(Maxrange-Xrange)/2])
end
set(gca,'ZLim', [ZLIM(1)-(Maxrange-Zrange)/2 ZLIM(2)+(Maxrange-Zrange)/2])
set(gca,'YLim', [YLIM(1)-(Maxrange-Yrange)/2 YLIM(2)+(Maxrange-Yrange)/2])
%set(gca,'Color',[0.2 0.2 0.2]);



%% Movie option
if PaperVersion2D==1
    pos=[-7600 0 -1300] - BaseShift;
else
    pos=[-5000 0 -2000] - BaseShift;
end

set(gca,'CameraPosition',pos);

if mov==1
    pause()
    pos_local=pos;
    for k = 1:40
        pos_local(1) = pos_local(1) + 100;
        pos_local(2) = pos_local(2) - 100;
        pos_local(3) = pos_local(3) + 100;
        set(gca,'CameraPosition',pos_local);
        %light('Position', pos_local, 'Style', 'infinite')
        pause(1)
        M(k) = getframe;
    end
end
%movie(M,5)


hold off
%colormap(MAP(1000,:))
%shading interp

%axis([-10 10 -10 10 -10 10])