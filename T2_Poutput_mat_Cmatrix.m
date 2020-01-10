clc
clear
save_directory='Beamforming��������';
mkdir(save_directory);
[fname,location]=uigetfile({'*.mat'},'mat�����ļ���ȡ','MultiSelect','off');%MultiSelect��ѡ
Sig = importdata(fullfile(location,fname)); %ѡ���ļ���������

%�����޶�Ƶ�ʹ۲췶Χ
df=Sig.freq(2)-Sig.freq(1);
maxFreq=round((Sig.rotor_speed/60*29)*2.5/df);
minFreq=round(1/3*(Sig.rotor_speed/60)/df);Freq=Sig.freq(minFreq):df:Sig.freq(maxFreq);

for k=minFreq:maxFreq  
cormatrix(:,:,k)=Sig.CC1(k,:).'*conj(Sig.CC1(k,:));
end


speed = 344;                              %����
zR = 0;                                 %�ؽ�����
zH = 0.4;                                 %��������

%% ������circle type������������������������������
% aperture size=0.2m,num of microphones is 360
theta=(0:length(Sig.fname)*12-1)*(360/length(Sig.fname)/12);
r=0.185;
xh=r*cos(theta*pi/180);yh=r*sin(theta*pi/180);Nxh=length(xh);
zh=zH*ones(1,length(xh));


%% ----------�ؽ���------------------------------------------
xr2=[-0.3:0.01:0.3];yr2=[-0.3:0.01:0.3];xrN=length(xr2);
xr=xr2'*ones(1,length(yr2));xr=reshape(xr,1,[]);
yr=yr2'*ones(1,length(xr2));yr=reshape(yr',1,[]);
zr=zR*ones(1,length(xr));Nxr=length(xr);
xrM=Nxr/xrN;    %--xrNΪ�еĳ��ȣ�xrMΪ�еĳ���


%% ---------�������þ�����ʽƵ�򷽷�����---------------

R=[xr' yr' zr'];
H=[xh' yh' zh'];
xrr=reshape(xr,xrM,xrN);
yrr=reshape(yr,xrM,xrN);
Dist=sqrt(R.^2*ones(size(H'))+ones(size(R))*(H').^2-2*R*H');
w=2*pi*Freq;


for k=1:length(w)
tic
%ww����ɨ������,cormatrix������ؾ���
%https://blog.csdn.net/qq_36300268/article/details/88739184 
%���˼·��������ȫ����һ��
A=exp(-1i*bsxfun(@times,Dist,w(k))/speed)/Nxh; %steer vector
pow(:,k)=sum(bsxfun(@times,conj(A),permute(cormatrix(:,:,k)'*A.',[2,1])),2);
pow1=reshape(pow(:,k),xrM,xrN);
h=figure
contourf(xrr,yrr,abs(pow1));
hold on
aplha=0:pi/40:2*pi;
r1=0.185;r2=0.13;
x1=r1*cos(aplha);x2=r2*cos(aplha);
y1=r1*sin(aplha);y2=r2*sin(aplha);
plot(x1,y1,'-');plot(x2,y2,'-');
axis equal
title(['Freq:',num2str(Freq(k)),'Hz'])
saveas(h,[save_directory,'\',strrep(char(fname),'.mat','-'),round(num2str(Freq(k))),'Hz','.png'])
close all
toc
end      

        
        




