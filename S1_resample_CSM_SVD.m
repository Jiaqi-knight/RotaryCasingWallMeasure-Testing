clc
clear
% close all
%% Parameter3
c = 344;                              %声速
zR = 0;                                 %重建距离
zH = 0.4;                                 %测量距离
%% 输出配置
save_directory='Beamforming基础版结果-改装版ver1';
mkdir(save_directory);
[fname,location]=uigetfile({'*.mat'},'mat参数文件读取','MultiSelect','off');%MultiSelect单选
Sig = importdata(fullfile(location,fname)); %选择文件导入数据
subfunction_path1='H:\Jiaqi-SJTU-DOIT\subfunction\subfunction_0_plot';
addpath(subfunction_path1);
% subfunction_path2='.\subfunction7';
% addpath(genpath(subfunction_path2));
%% 对Sig.DATA利用key信号进行锁相移动处理，为了使旋转的多段信号能够匹配
%新的代码思路放弃resample，采用随机cut掉一些点，保持每个都是最小
% Order_Cycle= 2048;   %每转采集点数目
tic
%figure
for k=1:length(Sig.DATA)
    [Pulse{k}]=keyRotation(Sig.DATA{1,k}(:,end),Sig.fs); 
    diff_Plus1(k)=min(diff(Pulse{k}));
end
ms=min(diff_Plus1);
for k=1:length(Sig.DATA)
    data_resample{1,k} = zeros(ms*(length(Pulse{k})-1),length(Sig.objectBF));
    for j = 1:length(Pulse{k})-1
        %对比两种处理办法：1.等间隔重采用 2.随机cut掉一些点 3. cut掉末尾的点
        data_resample{1,k}(ms*(j-1)+1:ms*j,:) = resample(Sig.DATA{1,k}(Pulse{k}(j):Pulse{k}(j+1)-1,:), ms, Pulse{k}(j+1)-Pulse{k}(j));
        data_resample_randcut{1,k}(ms*(j-1)+1:ms*j,:) = resample_randcut(Sig.DATA{1,k}(Pulse{k}(j):Pulse{k}(j+1)-1,:),ms);
        data_resample_endcut{1,k}(ms*(j-1)+1:ms*j,:) = resample_endcut(Sig.DATA{1,k}(Pulse{k}(j):Pulse{k}(j+1)-1,:),ms);
        
    end
%     plot(data_resample{1,k}(:,end-1));hold on;
    %对ref信号分析发现，与key信号完美契合，证明直接利用键向信号的可行性
end
toc
            L_signal =ms*(length(Pulse{k})-1);
            L_seg = round(L_signal/10);
            Wind = hamming(L_seg);
            Noverlap = round(L_seg/2);
            Nfft =2^(ceil(log2(L_seg))+1); 
for i_file=1:length(Sig.DATA)
            [the_freq12,data_freq12]=calculation_FFT(data_resample{i_file}(:,1:12),Sig.fs,1);
            [the_freq13,data_freq13]=calculation_FFT(data_resample{i_file}(:,1:12),Sig.fs,2);
            [the_freq22,data_freq22]=calculation_FFT(data_resample_randcut{i_file}(:,1:12),Sig.fs,1);
            [the_freq23,data_freq23]=calculation_FFT(data_resample_randcut{i_file}(:,1:12),Sig.fs,2);
            [the_freq32,data_freq32]=calculation_FFT(data_resample_endcut{i_file}(:,1:12),Sig.fs,1);
            [the_freq33,data_freq33]=calculation_FFT(data_resample_endcut{i_file}(:,1:12),Sig.fs,2);
            for k=1:12
            [data_freq0,the_freq0] = cpsd(Sig.DATA{i_file}(:,k),Sig.DATA{i_file}(:,13),Wind,Noverlap,Nfft,Sig.fs);
            CC0(:,i_file+length(Sig.DATA)*(k-1)) = data_freq0;      %CC1：通过ref消除
%             Data_reconstruct1(:,i_file+length(Sig.DATA)*(k-1)) = data_resample{i_file}(:,k);
%             Data_reconstruct2(:,i_file+length(Sig.DATA)*(k-1)) = data_resample_randcut{i_file}(:,k);
%             Data_reconstruct3(:,i_file+length(Sig.DATA)*(k-1)) = data_resample_endcut{i_file}(:,k);           
            tic
            [data_freq11,the_freq11] = cpsd(data_resample{i_file}(:,k),data_resample{i_file}(:,13),Wind,Noverlap,Nfft,Sig.fs);
            CC11(:,i_file+length(Sig.DATA)*(k-1)) = data_freq11;      %CC1：通过ref消除
            CC12(:,i_file+length(Sig.DATA)*(k-1)) =data_freq12(:,k);%已经进行key移相操作，是否CC1==CC2
            CC13(:,i_file+length(Sig.DATA)*(k-1)) =data_freq13(:,k);%已经进行key移相操作，是否CC1==CC2
            toc
            tic
            [data_freq21,the_freq21] = cpsd(data_resample_randcut{i_file}(:,k),data_resample_randcut{i_file}(:,13),Wind,Noverlap,Nfft,Sig.fs);
            CC21(:,i_file+length(Sig.DATA)*(k-1)) = data_freq21;      %CC1：通过ref消除
            CC22(:,i_file+length(Sig.DATA)*(k-1)) =data_freq22(:,k);%已经进行key移相操作，是否CC1==CC2
            CC23(:,i_file+length(Sig.DATA)*(k-1)) =data_freq23(:,k);%已经进行key移相操作，是否CC1==CC2
            toc
            tic
            [data_freq31,the_freq31] = cpsd(data_resample_endcut{i_file}(:,k),data_resample_endcut{i_file}(:,13),Wind,Noverlap,Nfft,Sig.fs);
            CC31(:,i_file+length(Sig.DATA)*(k-1)) = data_freq31;      %CC1：通过ref消除
            CC32(:,i_file+length(Sig.DATA)*(k-1)) =data_freq32(:,k);%已经进行key移相操作，是否CC1==CC2
            CC33(:,i_file+length(Sig.DATA)*(k-1)) =data_freq33(:,k);%已经进行key移相操作，是否CC1==CC2
            toc
            end
end

%对比两种处理办法：1.等间隔重采用 2.随机cut掉一些点 3. cut掉末尾的点
%1.cpsd-ref 2.fft-method1  3.fft-method2
%选择一个frequency
% frequency=3000;
% %index_freq0=round(frequency/(the_freq0(2)-the_freq0(1)))+1;
% %index_freq11=round(frequency/(the_freq11(2)-the_freq11(1)))+1;
% index_freq12=round(frequency/(the_freq12(2)-the_freq12(1)))+1;
% index_freq13=round(frequency/(the_freq13(2)-the_freq13(1)))+1;
% %index_freq21=round(frequency/(the_freq21(2)-the_freq21(1)))+1;
% index_freq22=round(frequency/(the_freq22(2)-the_freq22(1)))+1;
% index_freq23=round(frequency/(the_freq23(2)-the_freq23(1)))+1;
% %index_freq31=round(frequency/(the_freq31(2)-the_freq31(1)))+1;
% index_freq32=round(frequency/(the_freq32(2)-the_freq32(1)))+1;
% index_freq33=round(frequency/(the_freq33(2)-the_freq33(1)))+1;

% p=[];
% p=data_freq33(:,:)';
% % spectral matrix of the microphone measurements
% Spp = p*p'; % 30*30 
% Spp = (Spp + Spp')./2;
% 
% figure,
% subplot(211);
% imagesc(real(Spp)); % 
% title('spectral matrix of the measurements');
% subplot(212);
% stem(svd(Spp));
% 
% 
% 
% frequency=[2900;]
% for k=1:31
% R = CSM(Data_reconstruct', frequency,length(Data_reconstruct)/Sig.testPeriod,10035,10);
% figure
% stem(svd(R))
% end
CC{1}=CC0;CC{2}=CC11;CC{3}=CC12;CC{4}=CC13;CC{5}=CC21;CC{6}=CC22;CC{7}=CC23;CC{8}=CC31;CC{9}=CC32;CC{10}=CC33;
freq{1}=the_freq0';freq{2}=the_freq11';freq{3}=the_freq12;freq{4}=the_freq13;freq{5}=the_freq21';freq{6}=the_freq22;freq{7}=the_freq23;freq{8}=the_freq31';freq{9}=the_freq32;freq{10}=the_freq33;
for kc=1:10
a_mf=[];  
nk=size(CC{kc},2);
m=-nk/2:nk/2;
    for k=1:length(m)
        a_mf(:,k)=1/nk*CC{kc}(:,1:nk)*exp(m(k)*i*2*pi*(1:nk)/nk).'; 
    end
    GAMMA = 10*log10(abs(a_mf)/4e-10);
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    colormap(jet);
    imagesc(m,freq{kc},GAMMA); 
    axis xy; %ylim([1,30]);
end       
Freq1=[2000:500:20000];
for kf=1:length(Freq1)
%必须限定频率观察范围
h=figure;set(gcf,'Position',get(0,'ScreenSize'));%imshow();
matrix=[1,2,3,4,6,7,8,10,11,12]
matrixName={'Time:None+Freq:Cpsd';...
    'Time:resample+Freq:Cpsd';'Time:resample+Freq:fft(noSum)';'Time:resample+Freq:fft(Sum)';...
    'Time:randcut+Freq:Cpsd';'Time:randcut+Freq:fft(noSum)';'Time:randcut+Freq:fft(Sum)';...
    'Time:endcut+Freq:Cpsd';'Time:endcut+Freq:fft(noSum)';'Time:endcut+Freq:fft(Sum)';}

for kc=1:10
subplot(3,4,matrix(kc))   
df=freq{kc}(2)-freq{kc}(1);
oneFreq=round(Freq1(kf)/df);
Freq=freq{kc}(oneFreq)
cormatrix=[]
for k=1:length(Freq)
cormatrix(:,:,k)=CC{kc}(oneFreq,:).'*conj(CC{kc}(oneFreq,:));
end


speed = 344;                              %声速
zR = 0;                                 %重建距离
zH = 0.4;                                 %测量距离

%% ———circle type———————————————
% aperture size=0.2m,num of microphones is 360
theta=(0:length(Sig.fname)*12-1)*(360/length(Sig.fname)/12);
r=0.185;
xh=r*cos(theta*pi/180);yh=r*sin(theta*pi/180);Nxh=length(xh);
zh=zH*ones(1,length(xh));


%% ----------重建面------------------------------------------
xr2=[-0.3:0.01:0.3];yr2=[-0.3:0.01:0.3];xrN=length(xr2);
xr=xr2'*ones(1,length(yr2));xr=reshape(xr,1,[]);
yr=yr2'*ones(1,length(xr2));yr=reshape(yr',1,[]);
zr=zR*ones(1,length(xr));Nxr=length(xr);
xrM=Nxr/xrN;    %--xrN为列的长度，xrM为行的长度


%% ---------以下是用矩阵形式频域方法运算---------------

R=[xr' yr' zr'];
H=[xh' yh' zh'];
xrr=reshape(xr,xrM,xrN);
yrr=reshape(yr,xrM,xrN);
Dist=sqrt(R.^2*ones(size(H'))+ones(size(R))*(H').^2-2*R*H');
w=2*pi*Freq;


for k=1:length(w)
tic
%ww代表扫描向量,cormatrix代表相关矩阵
%https://blog.csdn.net/qq_36300268/article/details/88739184 
%借鉴思路，但不完全和他一致
A=exp(-1i*bsxfun(@times,Dist,w(k))/speed)/Nxh; %steer vector
pow(:,k)=sum(bsxfun(@times,conj(A),permute(cormatrix(:,:,k)'*A.',[2,1])),2);
pow1=reshape(pow(:,k),xrM,xrN);

contourf(xrr,yrr,abs(pow1));
hold on
aplha=0:pi/40:2*pi;
r1=0.185;r2=0.13;
x1=r1*cos(aplha);x2=r2*cos(aplha);
y1=r1*sin(aplha);y2=r2*sin(aplha);
plot(x1,y1,'-');plot(x2,y2,'-');
axis equal
title({['Freq:',num2str(Freq(k)),'Hz'];[matrixName{kc}]})
% saveas(h,[save_directory,'\',strrep(char(fname),'.mat','-'),round(num2str(Freq(k))),'Hz','.png'])
% close all
toc
end      



end
    
saveas(h,[save_directory,'\VS-',strrep(char(fname),'.mat','-'),round(num2str(Freq(k))),'Hz','.png'])
saveas(h,[save_directory,'\VS-',strrep(char(fname),'.mat','-'),round(num2str(Freq(k))),'Hz','.fig'])
end

       
function data=resample_randcut(data,ms)
data(randperm(length(data),length(data)-ms),:)=[];
if length(data)-ms>10
disp('速度均匀度太差');
end
end
function data1=resample_endcut(data,ms)
data1=data(1:ms,:);
if length(data)-ms>10
disp('速度均匀度太差');
end
end




