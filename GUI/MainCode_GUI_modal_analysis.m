function [] = GUI_modal()

SCR = get(0,'Screensize');  % Get screensize.
S.fh = figure('units','pixels','position',[SCR(1)+SCR(3)*0.025 ,SCR(2)+SCR(4)*0.055 , SCR(3)*0.95, SCR(4)*0.85],'name','GUI_1','resize','off');

S.frame1=uicontrol('style','frame','units','pixels','position',[450/1920*SCR(3) 700/1080*SCR(4) 870/1920*SCR(3) 150/1080*SCR(4)]);

S.txtotal=uicontrol('style','text','units','pixels','position',[470/1920*SCR(3) 710/1080*SCR(4) 60/1920*SCR(3) 130/1080*SCR(4)],'string','实验基础参数记录','value',0,'enable','on','fontsize',12,'fontweight','bold');
S.txtotal=uicontrol('style','text','units','pixels','position',[570/1920*SCR(3) 740/1080*SCR(4) 120/1920*SCR(3) 30/1080*SCR(4)],'string','实验通道数','value',0,'enable','on','fontsize',12,'fontweight','bold');
S.chtotal1 = uicontrol('style','edit','units','pixels','position',[540/1920*SCR(3) 710/1080*SCR(4) 150/1920*SCR(3) 30/1080*SCR(4)], 'string','输入整数', 'value',0, 'enable','on','fontsize',12,'foregroundcolor','red');

S.txtotal=uicontrol('style','text','units','pixels','position',[750/1920*SCR(3) 810/1080*SCR(4) 120/1920*SCR(3) 30/1080*SCR(4)],'string','采样频率','value',0,'enable','on','fontsize',12,'fontweight','bold');
S.chtotal2 = uicontrol('style','edit','units','pixels','position',[730/1920*SCR(3) 780/1080*SCR(4) 150/1920*SCR(3) 30/1080*SCR(4)], 'string','输入整数', 'value',0, 'enable','on','fontsize',12,'foregroundcolor','red');

S.txtotal=uicontrol('style','text','units','pixels','position',[750/1920*SCR(3) 740/1080*SCR(4) 130/1920*SCR(3) 30/1080*SCR(4)],'string','转子叶片数','value',0,'enable','on','fontsize',12,'fontweight','bold');
S.chtotal3 = uicontrol('style','edit','units','pixels','position',[730/1920*SCR(3) 710/1080*SCR(4) 150/1920*SCR(3) 30/1080*SCR(4)], 'string','输入实数', 'value',0, 'enable','on','fontsize',12,'foregroundcolor','red');

S.ch2 = uicontrol('style','edit','units','pixels','position',[550/1920*SCR(3) 780/1080*SCR(4) 150/1920*SCR(3) 30/1080*SCR(4)], 'string','输入有理数', 'value',0, 'enable','on','fontsize',12,'foregroundcolor','red');
S.tx2=uicontrol('style','text','units','pixels','position',[530/1920*SCR(3) 810/1080*SCR(4) 180/1920*SCR(3) 30/1080*SCR(4)],'string','阵列半径(mm)','value',0,'enable','on','fontweight','bold','fontsize',12,'fontweight','bold');
%% 界面设置模块1
S.ax2 = axes('units','pixels','position',[100/1920*SCR(3) 600/1080*SCR(4) 300/1920*SCR(3) 300/1080*SCR(4)]);
S.pb2 = uicontrol('style','pushbutton','units','pixels','position',[100/1920*SCR(3) 520/1080*SCR(4) 300/1920*SCR(3) 40/1080*SCR(4)],'string','Generate Transducer Array','fontsize',9,'fontweight','bold');
S.fi2 = uicontrol('style','toggle','units','pixels','position',[1250/1920*SCR(3) 805/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)], 'string','打开', 'value',0, 'enable','on','fontweight','bold');
S.finame2 = uicontrol('style','edit','units','pixels','position',[950/1920*SCR(3) 805/1080*SCR(4) 270/1920*SCR(3) 40/1080*SCR(4)],'string','传感器角度文件','value',0,'enable','on');
S.fi2b = uicontrol('style','toggle','units','pixels','position',[1250/1920*SCR(3) 755/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)], 'string','打开', 'value',0, 'enable','on','fontweight','bold');
S.finame2b = uicontrol('style','edit','units','pixels','position',[950/1920*SCR(3) 755/1080*SCR(4) 270/1920*SCR(3) 40/1080*SCR(4)],'string','传感器与采集通道次序文件','value',0,'enable','on');
S.fi = uicontrol('style','toggle','units','pixels','position',[1250/1920*SCR(3) 705/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)],'string','打开','value',0,'enable','on','fontweight','bold');
S.finame = uicontrol('style','edit','units','pixels','position',[950/1920*SCR(3) 705/1080*SCR(4) 270/1920*SCR(3) 40/1080*SCR(4)],'string','数据路径','value',0,'enable','on');
set(S.pb2,'callback',{@pb2_call,S})
set(S.fi2,'callback',{@fi2_call,S})
set(S.fi2b,'callback',{@fi2b_call,S})

%% 界面设置模块2
S.ax = axes('units','pixels','position',[70/1920*SCR(3) 160/1080*SCR(4) 400/1920*SCR(3) 300/1080*SCR(4)]);
S.pb1 = uicontrol('style','pushbutton', 'units','pixels','position',[100/1920*SCR(3) 5/1080*SCR(4) 300/1920*SCR(3) 40/1080*SCR(4)], 'string','Plot 1D Signal', 'fontsize',12,'fontweight','bold');
S.tg(1) = uicontrol('style','toggle','units','pixels','position',[120/1920*SCR(3) 50/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)],'string','TIME','value',0,'enable','off','fontweight','bold');
S.tg(2) = uicontrol('style','toggle','units','pixels','position',[185/1920*SCR(3) 50/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)],'string','FFT','value',0,'enable','off','fontweight','bold');
S.tg(3) = uicontrol('style','toggle','units','pixels','position',[250/1920*SCR(3) 50/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)],'string','STFT','value',0,'enable','off','fontweight','bold');
S.tg(4) = uicontrol('style','toggle','units','pixels','position',[315/1920*SCR(3) 50/1080*SCR(4) 60/1920*SCR(3) 40/1080*SCR(4)],'string','ABOUT','value',0,'enable','on','fontweight','bold');
S.tx = uicontrol('style','text','units','pixels','position',[70/1920*SCR(3) 160/1080*SCR(4) 360/1920*SCR(3) 300/1080*SCR(4)],'visible','off','string',{' ','This is a GUI with', 'DuctMode Vplot.','Hope you enjoy.',' ',' ','Copyright:','JiaqiWang 2019'},'fontsize',20,'fontweight','bold');
S.ch = uicontrol('style','edit','units','pixels','position',[440/1920*SCR(3) 465/1080*SCR(4) 30/1920*SCR(3) 30/1080*SCR(4)],'string','1','value',0,'enable','on');
set(S.pb1,'callback',{@pb_call,S})  % Set the callbacks.
set(S.fi,'callback',{@fi_call,S})  % Set the callbacks.
set(S.tg(:),{'callback'},{{@tg_call,S}})
%% 界面设置模块3
S.ax3 = axes('units','pixels','position',[580/1920*SCR(3) 120/1080*SCR(4) 640/1920*SCR(3) 520/1080*SCR(4)],'color','white');
S.ax4 = axes('units','pixels','position',[1400/1920*SCR(3) 540/1080*SCR(4) 400/1920*SCR(3) 300/1080*SCR(4)],'color','white');
S.ax5 = axes('units','pixels','position',[1400/1920*SCR(3) 50/1080*SCR(4) 400/1920*SCR(3) 300/1080*SCR(4)],'color','white');
% S.tx4 = uicontrol('style','text','units','pixels','position',[1450/1920*SCR(3) 900/1080*SCR(4) 300/1920*SCR(3) 20/1080*SCR(4)],'string','1 BPF','fontsize',9,'fontweight','bold');
% S.tx5 = uicontrol('style','text','units','pixels','position',[1450/1920*SCR(3) 400/1080*SCR(4) 300/1920*SCR(3) 20/1080*SCR(4)],'string','2 BPF','fontsize',9,'fontweight','bold');
% S.pb3b = uicontrol('style','pushbutton','units','pixels','position',[700/1920*SCR(3) 5/1080*SCR(4) 100/1920*SCR(3) 30/1080*SCR(4)],'string','Clear','fontsize',12,'fontweight','bold');
S.pb3 = uicontrol('style','pushbutton','units','pixels','position',[750/1920*SCR(3) 5/1080*SCR(4) 300/1920*SCR(3) 30/1080*SCR(4)],'string','Cut Off Modal Analysis','fontsize',12,'fontweight','bold');
set(S.pb3,'callback',{@pb3_call,S})

%% 界面1的函数
function [] = fi2_call(varargin)
global fname2 location2;
[fname2,location2]=uigetfile({'*.xlsx';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
set(varargin{3}.finame2,'string',fullfile(location2,fname2)) % Set edit to current slider.
function [] = fi2b_call(varargin)
global fname2b location2b;
[fname2b,location2b]=uigetfile({'*.xlsx';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
set(varargin{3}.finame2b,'string',fullfile(location2b,fname2b)) % Set edit to current slider.

function [] = pb2_call(varargin)
global fname2 location2 datum;
% Callback for pushbutton.
if strcmp(varargin{3}.pb2.String,'Generate Transducer Array')
    datum= xlsread(fullfile(location2,fname2));
    axes(varargin{3}.ax2);
    [trans_number] = transducer_array_graph(datum, str2num(get(varargin{3}.ch2,'string')));
end
%% 界面2的函数

function [] = fi_call(varargin)
global DATA Fs location2b fname2b fname location
Fs=str2num(varargin{3}.chtotal2.String);
% Callback for pushbutton.
[fname,location]=uigetfile({'*.wav';'*.*t'},'MultiSelect','on');%MultiSelect单选
set(varargin{3}.finame,'string',fullfile(location,fname)) % Set edit to current slider.
[data,Fs]= audioread(fullfile(location,fname),'native'); %选择文件导入数据
set(varargin{3}.tg,'enable','on');  % Turn on 'Fit' tab.
mic_transform_ab= xlsread(fullfile(location2b,fname2b));
DATA=(double(data(:,mic_transform_ab(:,2)))/(2^15)*10)/0.05;
save_directory = ['matData',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
DATA=(double(data(:,mic_transform_ab(:,2)))/(2^15)*10)/0.05;

function [] = pb_call(varargin)
global DATA Fs fname
Fs=str2num(varargin{3}.chtotal2.String);
axes(varargin{3}.ax);
if strcmp(varargin{3}.pb1.String,'Plot 1D Signal')
    yaxismax = max(max(abs(DATA))*1.2);
    plot(DATA(:,str2num(get(varargin{3}.ch,'string'))),'color',[0 0.4470 0.7410] );
    ylim([-yaxismax yaxismax])
    xlim([0,round(length(DATA)/Fs)*Fs]);
    xlimits = get(gca, 'Xlim'); % get the x-axis limits
    xlimits_time = xlimits./Fs;
    x2_time = ceil(xlimits_time(2));
    if x2_time>10;x2_time=11;end
    x3 = linspace(xlimits(1),xlimits(2),x2_time+1);
    x3_time = round(linspace(xlimits_time(1),xlimits_time(2),x2_time+1)./0.1)*0.1;
    set(gca,'xtick',x3)
    set(gca,'xticklabel',x3_time)
    xlabel('Time (s)','fontsize',12);
    ylabel('Amplitude','fontsize',12);
        title(fname,'fontsize',16);
elseif  strcmp(varargin{3}.pb1.String,'Plot 1D Spectrum')
    [the_freq,freq_dB]=frequencyDomainPlot_dB(DATA(:,str2num(varargin{3}.ch.String)),Fs,2.56);
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Amplitude (dB)','fontsize',12);
    title(fname,'fontsize',16);
elseif  strcmp(varargin{3}.pb1.String,'Plot 1D STFT')
    window =1024;% str2num(get(handles.edit11,'string'));
    noverlap =512;% str2num(get(handles.edit12,'string'));
    window = round(window);
    noverlap = round(noverlap);
    nfft = window;
    
    [S,F,T,P] =spectrogram(DATA(:,str2num(varargin{3}.ch.String)), window, noverlap, nfft, Fs);
    surf(T,F,10*log10(P),'edgecolor','none');
    axis tight;
    view(0,90);
    xlabel('Time (s)','fontsize',12);
    ylabel('Frequency (Hz)','fontsize',12);
        title(fname,'fontsize',16);
colorbar;
end

function [] = tg_call(varargin)
[h,S] = varargin{[1,3]};  % Get calling handle ans structure.
if get(h,'val')==0  % Here the Toggle is already pressed.
    set(h,'val',1) % To keep the Tab-like functioning.
end
L = get(S.ax,'children');  % The line object.
switch h
    case S.tg(1) %这是第一个界面(TIME)
        set(varargin{3}.pb1,'string','Plot 1D Signal')
        set(S.tg([2,3,4]),'val',0)
        set([S.ax,S.pb1,S.ch],{'visible'},{'on'})
        set([L],{'visible'},{'on'})
        set([S.tx;],{'visible'},{'off'})
        
    case S.tg(2)%这是第二个界面(FFT)
        set([S.ax,S.pb1,S.ch],{'visible'},{'on'})
        set(S.tg([1,3,4]),{'val'},{0})
        set([L],{'visible'},{'off'})
        set([S.tx],{'visible'},{'off'})
        set(varargin{3}.pb1,'string','Plot 1D Spectrum')
        
    case S.tg(3)%这是第三个界面(STFT)
        set([S.ax,S.pb1,S.ch],{'visible'},{'on'})
        set(S.tg([1,2,4]),{'val'},{0})
        set([L],{'visible'},{'off'})
        set([S.tx],{'visible'},{'off'})
        set(varargin{3}.pb1,'string','Plot 1D STFT')
        
    otherwise
        set(S.tg([1,2,3]),{'val'},{0})
        set(S.tx,'visible','on')
        set([S.ax;S.pb1;L;S.ch],{'visible'},{'off'})
end
%% 界面3的函数

function [] = pb3_call(varargin)
global fname location nk Fs Blade datum location2b fname2b;

Fs=str2num(varargin{3}.chtotal2.String);
Blade=str2num(varargin{3}.chtotal3.String);
[data,Fs]= audioread(fullfile(location,fname),'native'); %选择文件导入数据
mic_transform_ab=xlsread(fullfile(location2b,fname2b));
save_directory = ['matData',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
DATA=(double(data(:,mic_transform_ab(:,2)))/(2^15)*10)/0.05;

if strcmp(varargin{3}.pb3.String,'Cut Off Modal Analysis')
    Phase=datum/180*pi;
        DATA=  DATA./rms(  DATA);%对信号做正则化处理,正则化稀疏度更好
    [the_freq1,freq_dB1]=frequencyDomainPlot_dB2( DATA,Fs,2.56);
    nk=str2num(varargin{3}.chtotal1.String);
    BPF=the_freq1(find(freq_dB1(:,1)==max(freq_dB1(:,1))));rotor_speed=BPF/Blade;
    listName={'1BPF';'2BPF';};
    axes(varargin{3}.ax3);
        cla(varargin{3}.ax3,'reset');
    [GAMMA,freq]=wavemode_calculation_CPSD_phase(  DATA,Fs,nk,Phase,rotor_speed*60,save_directory,1,fname);
    df=freq(2)-freq(1);list=[round((Blade-0.5)*rotor_speed/df) round((Blade+0.5)*rotor_speed/df);round((2*Blade-0.5)*rotor_speed/df) round((2*Blade+0.5)*rotor_speed/df);];
    for k=1:2
        [x(k),y(k)]=find(GAMMA==max(max(GAMMA(list(k,1):list(k,2),:))),1);Zdata=GAMMA(x(k),:);Zdata(y(k))=0;[y_1(k)]=find(Zdata==max(Zdata),1);
        o_RI(k)=x(k)*df;m_RI(k)=y(k)-(nk/2+1);m_RI_1(k)=y_1(k)-(nk/2+1);
        text(m_RI(k),o_RI(k),{[num2str(round(m_RI(k)*10)/10),',',num2str(round(GAMMA(x(k),y(k))))];[num2str(o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60))]});text(m_RI_1(k),o_RI(k),[num2str(m_RI_1(k)),',',num2str(round(GAMMA(x(k),y_1(k))))]);
    end
    for k=1:2
        if k==1
            axes(varargin{3}.ax4);
                    cla(varargin{3}.ax4,'reset');
        else
            axes(varargin{3}.ax5);
                    cla(varargin{3}.ax5,'reset');
        end
        bar([-nk/2:nk/2],GAMMA(x(k),:));hold on;
        ylim([40 70]);
        title([char(fname),'-',num2str(round(rotor_speed*60)),'rpm','-',listName{k}],'FontSize',14);
            xlabel('Mode Number：m','FontSize',12);ylabel('Amplitude (dB)','FontSize',12);
    end
end

