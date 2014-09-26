%% Script for example confusion matrix
L = {'Danger-F' 'Danger-M' 'Love-F' 'Love-M' 'Hunger-F' 'Hunger-M'};
CLim=[0 0.5];

%create a Non selective invariant cell
p=1/12;
NS_Inv=[p p 0 0 0 0;p p 0 0 0 0;0 0 p p 0 0;0 0 p p 0 0;0 0 0 0 p p;0 0 0 0 p p];

figure1 = figure('PaperSize',[20.98404194812 29.67743169791],...
     'Colormap',[0.0416666679084301 0 0;0.201388895511627 0 0;0.361111104488373 0 0;0.520833313465118 0 0;0.680555582046509 0 0;0.840277791023254 0 0;1 0 0;1 0.125 0;1 0.25 0;1 0.375 0;1 0.5 0;1 0.625 0;1 0.75 0;1 0.875 0;1 1 0;1 1 0.0232493001967669;1 1 0.0464986003935337;1 1 0.0697479024529457;1 1 0.0929972007870674;1 1 0.11624650657177;1 1 0.139495804905891;1 1 0.209243699908257;1 1 0.278991609811783;1 1 0.348739504814148;1 1 0.418487399816513;1 1 0.488235294818878;1 1 0.557983219623566;1 1 0.627731084823608;1 1 0.697479009628296;1 1 0.767226874828339;1 1 0.836974799633026;1 1 0.906722664833069;1 1 0.976470589637756;1 1 0.977229595184326;1 1 0.977988600730896;1 1 0.978747606277466;1 1 0.97950667142868;1 1 0.98026567697525;1 1 0.98102468252182;1 1 0.98178368806839;1 1 0.98254269361496;1 1 0.98330169916153;1 1 0.984060704708099;1 1 0.984819710254669;1 1 0.985578775405884;1 1 0.986337780952454;1 1 0.987096786499023;1 1 0.987855792045593;1 1 0.988614797592163;1 1 0.989373803138733;1 1 0.990132808685303;1 1 0.990891814231873;1 1 0.991650879383087;1 1 0.992409884929657;1 1 0.993168890476227;1 1 0.993927896022797;1 1 0.994686901569366;1 1 0.995445907115936;1 1 0.996204912662506;1 1 0.996963918209076;1 1 0.997722983360291;1 1 0.99848198890686;1 1 0.99924099445343;1 1 1]);

 % Create axes
Ticks = 1:length(L);
axes1 = axes('Parent',figure1,'YTick', Ticks, 'YTickLabel',L,...
    'YDir','reverse',...
    'XTick', Ticks, 'XTickLabel',L,...
    'Layer','top',...
    'CLim',CLim, 'TickLength',[0 0]);
hold(axes1,'all');
% Create image
image(NS_Inv,'Parent',axes1,'CDataMapping','scaled');
axis([0.5 length(L)+0.5 0.5 length(L)+0.5])

% Create colorbar
colorbar('peer',axes1);

%add some legend
xlabel('Predicted Vocalization');
ylabel('Actual Vocalization');

%create a Non selective non-invariant cell
p=1/6;
NS_NInv=[p 0 0 0 0 0;0 p 0 0 0 0;0 0 p 0 0 0;0 0 0 p 0 0;0 0 0 0 p 0;0 0 0 0 0 p];

figure2 = figure('PaperSize',[20.98404194812 29.67743169791],...
     'Colormap',[0.0416666679084301 0 0;0.201388895511627 0 0;0.361111104488373 0 0;0.520833313465118 0 0;0.680555582046509 0 0;0.840277791023254 0 0;1 0 0;1 0.125 0;1 0.25 0;1 0.375 0;1 0.5 0;1 0.625 0;1 0.75 0;1 0.875 0;1 1 0;1 1 0.0232493001967669;1 1 0.0464986003935337;1 1 0.0697479024529457;1 1 0.0929972007870674;1 1 0.11624650657177;1 1 0.139495804905891;1 1 0.209243699908257;1 1 0.278991609811783;1 1 0.348739504814148;1 1 0.418487399816513;1 1 0.488235294818878;1 1 0.557983219623566;1 1 0.627731084823608;1 1 0.697479009628296;1 1 0.767226874828339;1 1 0.836974799633026;1 1 0.906722664833069;1 1 0.976470589637756;1 1 0.977229595184326;1 1 0.977988600730896;1 1 0.978747606277466;1 1 0.97950667142868;1 1 0.98026567697525;1 1 0.98102468252182;1 1 0.98178368806839;1 1 0.98254269361496;1 1 0.98330169916153;1 1 0.984060704708099;1 1 0.984819710254669;1 1 0.985578775405884;1 1 0.986337780952454;1 1 0.987096786499023;1 1 0.987855792045593;1 1 0.988614797592163;1 1 0.989373803138733;1 1 0.990132808685303;1 1 0.990891814231873;1 1 0.991650879383087;1 1 0.992409884929657;1 1 0.993168890476227;1 1 0.993927896022797;1 1 0.994686901569366;1 1 0.995445907115936;1 1 0.996204912662506;1 1 0.996963918209076;1 1 0.997722983360291;1 1 0.99848198890686;1 1 0.99924099445343;1 1 1]);

 % Create axes
Ticks = 1:length(L);
axes2 = axes('Parent',figure2,'YTick', Ticks, 'YTickLabel',L,...
    'YDir','reverse',...
    'XTick', Ticks, 'XTickLabel',L,...
    'Layer','top',...
    'CLim',CLim, 'TickLength',[0 0]);
hold(axes2,'all');
% Create image
image(NS_NInv,'Parent',axes2,'CDataMapping','scaled');
axis([0.5 length(L)+0.5 0.5 length(L)+0.5])

% Create colorbar
colorbar('peer',axes2);

%add some legend
xlabel('Predicted Vocalization');
ylabel('Actual Vocalization');


%create a Selective non-invariant cell
p=1/6;
p2=1/36;
S_NInv=[p2 p2 p2 p2 p2 p2;p2 p2 p2 p2 p2 p2;0 0 p 0 0 0;0 0 0 p 0 0;p2 p2 p2 p2 p2 p2;p2 p2 p2 p2 p2 p2];

figure3 = figure('PaperSize',[20.98404194812 29.67743169791],...
     'Colormap',[0.0416666679084301 0 0;0.201388895511627 0 0;0.361111104488373 0 0;0.520833313465118 0 0;0.680555582046509 0 0;0.840277791023254 0 0;1 0 0;1 0.125 0;1 0.25 0;1 0.375 0;1 0.5 0;1 0.625 0;1 0.75 0;1 0.875 0;1 1 0;1 1 0.0232493001967669;1 1 0.0464986003935337;1 1 0.0697479024529457;1 1 0.0929972007870674;1 1 0.11624650657177;1 1 0.139495804905891;1 1 0.209243699908257;1 1 0.278991609811783;1 1 0.348739504814148;1 1 0.418487399816513;1 1 0.488235294818878;1 1 0.557983219623566;1 1 0.627731084823608;1 1 0.697479009628296;1 1 0.767226874828339;1 1 0.836974799633026;1 1 0.906722664833069;1 1 0.976470589637756;1 1 0.977229595184326;1 1 0.977988600730896;1 1 0.978747606277466;1 1 0.97950667142868;1 1 0.98026567697525;1 1 0.98102468252182;1 1 0.98178368806839;1 1 0.98254269361496;1 1 0.98330169916153;1 1 0.984060704708099;1 1 0.984819710254669;1 1 0.985578775405884;1 1 0.986337780952454;1 1 0.987096786499023;1 1 0.987855792045593;1 1 0.988614797592163;1 1 0.989373803138733;1 1 0.990132808685303;1 1 0.990891814231873;1 1 0.991650879383087;1 1 0.992409884929657;1 1 0.993168890476227;1 1 0.993927896022797;1 1 0.994686901569366;1 1 0.995445907115936;1 1 0.996204912662506;1 1 0.996963918209076;1 1 0.997722983360291;1 1 0.99848198890686;1 1 0.99924099445343;1 1 1]);

 % Create axes
Ticks = 1:length(L);
axes3 = axes('Parent',figure3,'YTick', Ticks, 'YTickLabel',L,...
    'YDir','reverse',...
    'XTick', Ticks, 'XTickLabel',L,...
    'Layer','top',...
    'CLim',CLim, 'TickLength',[0 0]);
hold(axes3,'all');
% Create image
image(S_NInv,'Parent',axes3,'CDataMapping','scaled');
axis([0.5 length(L)+0.5 0.5 length(L)+0.5])

% Create colorbar
colorbar('peer',axes3);

%add some legend
xlabel('Predicted Vocalization');
ylabel('Actual Vocalization');


%create a Selective invariant cell
p=1/12;
p2=1/36;
S_NInv=[p2 p2 p2 p2 p2 p2;p2 p2 p2 p2 p2 p2;0 0 p p 0 0;0 0 p p 0 0;p2 p2 p2 p2 p2 p2;p2 p2 p2 p2 p2 p2];

figure4 = figure('PaperSize',[20.98404194812 29.67743169791],...
     'Colormap',[0.0416666679084301 0 0;0.201388895511627 0 0;0.361111104488373 0 0;0.520833313465118 0 0;0.680555582046509 0 0;0.840277791023254 0 0;1 0 0;1 0.125 0;1 0.25 0;1 0.375 0;1 0.5 0;1 0.625 0;1 0.75 0;1 0.875 0;1 1 0;1 1 0.0232493001967669;1 1 0.0464986003935337;1 1 0.0697479024529457;1 1 0.0929972007870674;1 1 0.11624650657177;1 1 0.139495804905891;1 1 0.209243699908257;1 1 0.278991609811783;1 1 0.348739504814148;1 1 0.418487399816513;1 1 0.488235294818878;1 1 0.557983219623566;1 1 0.627731084823608;1 1 0.697479009628296;1 1 0.767226874828339;1 1 0.836974799633026;1 1 0.906722664833069;1 1 0.976470589637756;1 1 0.977229595184326;1 1 0.977988600730896;1 1 0.978747606277466;1 1 0.97950667142868;1 1 0.98026567697525;1 1 0.98102468252182;1 1 0.98178368806839;1 1 0.98254269361496;1 1 0.98330169916153;1 1 0.984060704708099;1 1 0.984819710254669;1 1 0.985578775405884;1 1 0.986337780952454;1 1 0.987096786499023;1 1 0.987855792045593;1 1 0.988614797592163;1 1 0.989373803138733;1 1 0.990132808685303;1 1 0.990891814231873;1 1 0.991650879383087;1 1 0.992409884929657;1 1 0.993168890476227;1 1 0.993927896022797;1 1 0.994686901569366;1 1 0.995445907115936;1 1 0.996204912662506;1 1 0.996963918209076;1 1 0.997722983360291;1 1 0.99848198890686;1 1 0.99924099445343;1 1 1]);

 % Create axes
Ticks = 1:length(L);
axes4 = axes('Parent',figure4,'YTick', Ticks, 'YTickLabel',L,...
    'YDir','reverse',...
    'XTick', Ticks, 'XTickLabel',L,...
    'Layer','top',...
    'CLim',CLim, 'TickLength',[0 0]);
hold(axes4,'all');
% Create image
image(S_NInv,'Parent',axes4,'CDataMapping','scaled');
axis([0.5 length(L)+0.5 0.5 length(L)+0.5])

% Create colorbar
colorbar('peer',axes4);

%add some legend
xlabel('Predicted Vocalization');
ylabel('Actual Vocalization');