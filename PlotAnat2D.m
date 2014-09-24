function []=PlotAnat2D(LR,RC,DV,VAR,Size,Color,Colorcat,Colorcatname,FIXED_max,x,max_dot)

%LR: Left-Right position
%RC: Rostrocaudal position
%DV: Dorsoventral position
%FIXED_max = 1;

%RC=-RC;

%% Getting VAR ready for plotting.
% If VAR is always a non-null
% integer then you can use it directly in cubehelix by specifying x=1
% otherwise, the function first multiplies VARY by 10^x to make sure its
% variability is maximized after ceiling it to an integer
if nargin<10
    x=3;% here I'm doing *1000 as an example
end
if nargin<11
    max_dot=16;% here I'm doing *1000 as an example
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

%% uper bound the diameters of the spheres to 8
VAR_diam = VAR_cube./(max(VAR_cube)./max_dot);


%% Plot
figure()
NU = length(LR); % this is the number of points you have in your plot
if Color && Size
    if isempty(Colorcat)
        for ss=1:NU
            plot3(LR(ss),RC(ss), DV(ss),'ko','MarkerEdgeColor','none','LineStyle','none','MarkerFaceColor',GRAD(VAR_cube(ss),:), 'MarkerSize',VAR_diam(ss));
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
            %GRAD2=[0.8 0.8 0.8; 0 1 1; 0 0.8 1; 1 0.3 0; 1 0 0; 1 0 0 ; 1 0 0.6; 1 0.2 0.5; 0 1 0; 0 0 0];
            GRAD2=[0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0; 0.8 0.8 0.8];
        elseif length(Colorcatname)==2
            GRAD2=[0.8 0.8 0.8; 0.7 0 0.7];
        elseif length(Colorcatname)==4
            GRAD2=[0 0.5 1;0 1 1;1 0.75 0;0.5 0 0.1];
        else
            fprintf(1, 'Problem with the GRAD2 definition, the length of the color code is %d and the colormap (GRAD2) can only be 9 or 10', length(Colorcatname)); 
        end
        for ss=1:NU
            plot3(LR(ss),RC(ss), DV(ss),'ko','MarkerEdgeColor','none','LineStyle','none','MarkerFaceColor',GRAD2(Colorcat(ss),:), 'MarkerSize',VAR_diam(ss));
            hold on
        end
        for ii=1:length(Colorcatname)
            text(-1000,1400,-100.*ii,Colorcatname(ii), 'Color', GRAD2(ii,:))
        end
    end
end
if Color && ~Size
    if isempty(Colorcat)
        for ss=1:NU
            plot3(LR(ss),RC(ss), DV(ss),'ko','MarkerEdgeColor','none','LineStyle','none','MarkerFaceColor',GRAD(VAR_cube(ss),:), 'MarkerSize',max(VAR_diam));
            hold on
        end
    else
        if length(Colorcatname)==9
            GRAD2=[0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
        elseif length(Colorcatname)==7
            %GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.8 0.8 0.8];%Great with dark font
            GRAD2=[1 1 0; 1 0 0; 1 0.6 0; 0 0.5 1; 0 1 1; 0 1 0 ; 0.3 0.3 0.3];%Great with white font
        elseif length(Colorcatname)==10
            GRAD2=[0.8 0.8 0.8; 0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
        elseif length(Colorcatname)==2
            %GRAD2=[0.8 0.8 0.8; 0.7 0 0.7];%Great with dark font
            GRAD2=[0.6 0.6 0.6; 0.7 0 0.7];%Great with white font
        elseif length(Colorcatname)==4
            GRAD2=[0 0.5 1;0 1 1;1 0.75 0;0.5 0 0.1];
        else
            fprintf(1, 'Problem with the GRAD2 definition, the length of the color code is %d and the colormap (GRAD2) can only be 9 or 10', length(Colorcatname)); 
        end
        for ss=1:NU
            plot3(LR(ss),RC(ss), DV(ss),'ko','MarkerEdgeColor','none','LineStyle','none','MarkerFaceColor',GRAD2(Colorcat(ss),:));%, 'MarkerSize',max(VAR_diam));
            hold on
        end
        for ii=1:length(Colorcatname)
            text(-1000,1400,-100.*ii,Colorcatname(ii), 'Color', GRAD2(ii,:))
        end
    end
       
end

if ~Color && Size
    for ss=1:NU
        plot3(LR(ss),RC(ss), DV(ss),'ko','MarkerEdgeColor','k','LineStyle','none','MarkerFaceColor','y', 'MarkerSize',VAR_diam(ss));
        hold on
    end
end
%set(gca,'Color',[0.2 0.2 0.2]);
pos=[-80000 0 0];    
set(gca,'CameraPosition',pos);
end
