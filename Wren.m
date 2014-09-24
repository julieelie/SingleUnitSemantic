function [data,stimProtocol]= Wren(block_name,tank_dir,parameters)
% FULLCOH

% INPUT
% block_name - string identifying block to analyze
% tank_dir - string identifying path of block to analyze
% parameters - cell of strings pertaining to stimuli parameters to parse
%   out

% OUTPUT
% data
%  .block = block name
%  .tank = tank name
%  .Fs = sampling frequency of data (defined at beginning of code)
%  .stimProtocol = cells containing vector of each column of
%    stimProtocol.txt file
%  .LFPdata = matrix where each row is a vector of LFP signal for a single
%    electrode
%  .stim{#} - cell where {#} pertains to stimuli of same number
%    .start_time = vector containing start time of each stimuli (s)
%    .end_time = vector containing end time of each stimuli (s)
%    .bg_start = vector containing duration from mid_point between previous 
%        stimuli and current stimuli (s)
%    .bg_end = vector containing duration from mid_point between next 
%        stimuli and current stimuli (s)
%    .trial_num = vector containing identifying trial ID #
%    .LFPrms = (meant to contain rms of each trial)
%  .stim_dict = stimuli identification numbers found in block (not all
%    stimuli occur in each block)
% specific_stims = vector of stimuli fitting parameters


Fs=381.4697;

%% initialize
if nargin==0
    [file_name,path_name]=uigetfile('.f32','Choose first LFP file','/auto/fdata/');
    name=[path_name,file_name];
    slashes=strfind(name,'/');
    if isempty(slashes)
        slashes=strfind(name,'\');
    end
    tank_dir=name(1:slashes(end)-1);
    spaces=strfind(name,' ');
    block_name=name(slashes(end)+1:spaces(1)-1);
elseif nargin<2
    tank_dir=uigetdir('/auto/fdata/','Choose tank directory');
end

%% record identifying information
slashes=strfind(tank_dir,'/');
cohdata.tank=tank_dir(slashes(end)+1:end);
disp(sprintf('tank recorded: ''%s''',cohdata.tank));
cohdata.block=block_name;
disp(sprintf('block recorded: ''%s''',cohdata.block));

%% load in LFP & epoch data
disp('loading in LFP data:');
if ~exist(sprintf('%s %s file.mat',cohdata.tank,cohdata.block),'file')
    disp(sprintf('Previously parsed file ''%s %s file.mat'' not found. Parsing original data...',cohdata.tank,cohdata.block));
    [data,epoch_data]=LFPparse(block_name,tank_dir,Fs);
else
    disp(sprintf('loading ''%s %s file.mat''...',cohdata.tank,cohdata.block))
    data=load(sprintf('%s %s file.mat',cohdata.tank,cohdata.block));
    % Load in epoch data
    if exist(sprintf('%s/%s epoc Stim.txt',tank_dir,block_name),'file')
        epoch_data=dlmread(sprintf('%s/%s epoc Stim.txt',tank_dir,block_name),'\t');
    elseif exist(sprintf('%s/%s epoc Stm+.txt',tank_dir,block_name),'file')
        epoch_data=dlmread(sprintf('%s/%s epoc Stm+.txt',tank_dir,block_name),'\t');
    else
        disp('Cannot find epoch file, exiting...');
    end
end


%% determine stimuli to analyze
if exist('parameters','var')
    if ~isempty(parameters)
        [specific_stims,cohdata.parameters] = specifyStim(data.stimProtocol,parameters);
    else
        specific_stims=data.stim_dict;
        cohdata.parameters='all';
    end
else
    [specific_stims,cohdata.parameters] = specifyStim(data.stimProtocol);
end
% remove stims that aren't in block
good=zeros(size(specific_stims));
for index=1:length(specific_stims)
    if find(specific_stims(index)==data.stim_dict)
        good(index)=1;
    end
end
specific_stims=specific_stims(logical(good));


end


%%
function [data,epoch_data]=LFPparse(block_name,tank_dir,Fs)
% LFPPARSE

save_data_file=1;

% Locate LFP files
if nargin==0
    [file_name,path_name]=uigetfile('.f32','Choose first LFP file','/auto/fdata/');
    name=[path_name,file_name];
    slashes=strfind(name,'/');
    if isempty(slashes)
        slashes=strfind(name,'\');
    end
    tank_dir=name(1:slashes(end)-1);
    spaces=strfind(name,' ');
    block_name=name(slashes(end)+1:spaces(1)-1);
    Fs=inputdlg('What is the sampling frequency? (Hz)','Fs',1,{'381.4697'});
    Fs=str2num(Fs{1});
elseif nargin<2
    tank_dir=uigetdir('/auto/fdata/','Choose tank directory');
    Fs=inputdlg('What is the sampling frequency? (Hz)','Fs',1,{'381.4697'});
    Fs=str2num(Fs{1});
end


% Record identifying info
data.block=block_name;
slashes=strfind(tank_dir,'/');
data.tank=tank_dir(slashes(end)+1:end);
data.Fs=Fs;


% Load in epoch data
if exist(sprintf('%s/%s epoc Stim.txt',tank_dir,block_name),'file')
    epoch_data=dlmread(sprintf('%s/%s epoc Stim.txt',tank_dir,block_name),'\t');
elseif exist(sprintf('%s/%s epoc Stm+.txt',tank_dir,block_name),'file')
    epoch_data=dlmread(sprintf('%s/%s epoc Stm+.txt',tank_dir,block_name),'\t');
else
    disp('Cannot find epoch file, exiting...');
end


% Load in stimuli dictionary
if exist(sprintf('%sStims/%s/stimProtocol.txt',tank_dir(1:slashes(end-1)),data.tank))
    fid=fopen(sprintf('%sStims/%s/stimProtocol.txt',tank_dir(1:slashes(end-1)),data.tank));
else
    [file_name,path_name]=uigetfile('.txt','Choose stimProtocol file','/auto/fdata/');
    fid=fopen([path_name file_name]);
end
data.stimProtocol=textscan(fid,'%f %s %s %s %s %s %s %s %s');
fclose(fid);


% Determine which electrodes' data exist
electrode_dict=[];
LFPfiles=ls(sprintf('%s/%s LFPs*',tank_dir,block_name));
index=strfind(LFPfiles,'LFPs');
total_electrodes=length(index);
for electrode_num=1:total_electrodes
    if LFPfiles(index+6)==' '
        electrode_dict=[electrode_dict;str2num(LFPfiles(index(electrode_num)+5))];
    else %two digit
        electrode_dict=[electrode_dict;str2num(LFPfiles(index(electrode_num)+5:index(electrode_num)+6))];
    end
end
electrode_dict=sort(electrode_dict);


% Load in LFP data for all electrodes
data.LFPdata=[];
for electrode_num=1:max(electrode_dict)
    file=fopen(sprintf('%s/%s LFPs %d stream.f32',tank_dir,block_name,electrode_dict(electrode_num)),'r');
    data.LFPdata(electrode_dict(electrode_num),:)=fread(file, inf, 'float32')';
    fclose(file);
end


%Process epoch data
data.stim=trial_data();
data.stim_dict=[];
for index=1:size(epoch_data,1)
    if index==1 %first trial
        prev_end=0;
        latter_start=epoch_data(index+1,2);
    elseif index==size(epoch_data,1) %last trial
        latter_start=epoch_data(end,3);
        prev_end=epoch_data(index-1,3);
    else
        prev_end=epoch_data(index-1,3);
        latter_start=epoch_data(index+1,2);
    end
    if nnz(data.stim_dict==epoch_data(index,1)) %not first time of this stimulus
        data.stim(epoch_data(index,1)).start_time=[data.stim(epoch_data(index,1)).start_time;epoch_data(index,2)];
        data.stim(epoch_data(index,1)).end_time=[data.stim(epoch_data(index,1)).end_time;epoch_data(index,3)];
        data.stim(epoch_data(index,1)).bg_start=[data.stim(epoch_data(index,1)).bg_start;(epoch_data(index,2)-prev_end)/2];
        data.stim(epoch_data(index,1)).bg_end=[data.stim(epoch_data(index,1)).bg_end;(latter_start-epoch_data(index,3))/2];
        data.stim(epoch_data(index,1)).trial_num=[data.stim(epoch_data(index,1)).trial_num;index];
        trial_num=length(data.stim(epoch_data(index,1)).trial_num);
        start_time=epoch_data(index,2)-(epoch_data(index,2)-prev_end)/2;
        end_time=epoch_data(index,3)+(latter_start-epoch_data(index,3))/2;

        data.stim(epoch_data(index,1)).LFPrms{trial_num}=std(data.LFPdata(:,start_time*Fs:end_time*Fs));

    else %first time of this stimulus
        data.stim(epoch_data(index,1)).start_time=epoch_data(index,2);
        data.stim(epoch_data(index,1)).end_time=epoch_data(index,3);
        data.stim(epoch_data(index,1)).bg_start=(epoch_data(index,2)-prev_end)/2;
        data.stim(epoch_data(index,1)).bg_end=(latter_start-epoch_data(index,3))/2;
        data.stim(epoch_data(index,1)).trial_num=index;
        trial_num=length(data.stim(epoch_data(index,1)).trial_num);
        start_time=epoch_data(index,2)-(epoch_data(index,2)-prev_end)/2;
        end_time=epoch_data(index,3)+(latter_start-epoch_data(index,3))/2;
        
        data.stim(epoch_data(index,1)).LFPrms{trial_num}=std(data.LFPdata(:,start_time*Fs:end_time*Fs));

        data.stim_dict=[data.stim_dict;epoch_data(index,1)];
    end

end

data.stim_dict=sort(data.stim_dict);


% Save data cell structure
if save_data_file
    save(sprintf('%s %s file.mat',data.tank,data.block),'-struct','data');
end

end


%%
function [specific_stims,input] = specifyStim(stimProtocol,input)
% SPECIFYSTIM - identifies specific stimuli based on input criteria
% 
% Input data: 
%   stimProtocol - loaded by LFPparse.m from specific stimProtocol.txt file
% Input parameters: a cell array of strings


%initialize values
indexfinal=ones(size(stimProtocol{1},1),1);
parameters=[];

%build list of specification parameters in stimprotocol file
values=[];
for type=4:9
    values = [values;unique(stimProtocol{type});'--'];
end

%toggle ui tool if no 'input' variable
if ~exist('input','var')
    index = listdlg('PromptString','Specify Parameters:','SelectionMode','multiple','ListString',values);
    input=values(index);
    input(strcmp(input,'--'))=[]; %remove any partitions if input
elseif ~isempty(input)
    if ~iscellstr(input)
    error('input parameter(s) needs to be a cell array of strings');
    return
    end
end

%do not perform discrimination of stimuli if input is an empty variable
if ~isempty(input)
    %match input discrimination specifications to parameter types
    discrim=[];
    for in=1:length(input)
        match=0;
        for type=4:9
            index=find(strcmp(input(in),unique(stimProtocol{type})));
            if index
                discrim=[discrim,type];
                match=match+1;
            end
        end
        if ~match
            error(sprintf('%s did not match any values in stimprotocol file',input(in)));
            return
        elseif match>1
            error(sprintf('%s matched more than one value in stimprotocol file',input(in)));
            return
        end
    end
    groups=unique(discrim);
    
    %discriminate stimuli for each parameter type
    for type=1:length(groups)
        index=zeros(size(stimProtocol{1},1),1);
        inputcurrent=find(discrim==groups(type));
        for in=1:length(inputcurrent)
            indexcurrent = strcmp(input(inputcurrent(in)),stimProtocol{groups(type)});
            index(indexcurrent)=1;
        end
        indexfinal(indexfinal==1 & index==1)=1;
        indexfinal(indexfinal==0 | index==0)=0;
    end
end

% determine specific stimuli to analyze
specific_stims=stimProtocol{1}(logical(indexfinal));

end


