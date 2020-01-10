%wjq-2019-11-18
%代码目的：前处理，准备数据，并存储自相关，互相关，DATA结果
clc
clear
% subfunction_path1='H:\动画预测一体化程序通用版_旋转机匣_ver1\subfunction\subfunction_1';
% addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
save_directory = ['波束成形结果',date];  %频谱图存储文件夹

if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
    L_signal =fs*testPeriod;
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft =2^(ceil(log2(L_seg))+1); 
tic
for i_file=1:length(fname)
    close all
    Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    DATA{i_file}=Data(:,objectBF);
    
        for k=1:12
            [temp,freq] = cpsd(DATA{i_file}(:,k),DATA{i_file}(:,13),Wind,Noverlap,Nfft,fs);
            [temp1,freq1] = cpsd(DATA{i_file}(:,1:12),DATA{i_file}(:,1:12),Wind,Noverlap,Nfft,fs);
            CC1(:,i_file+length(fname)*(k-1)) = temp;
            CC2{i_file} =temp1;
        end
end

[Pulse,rotor_speed]=keyRotation(DATA{1}(:,end),fs);
h=figure     
plot(freq1,10*log10(sum(cell2mat(CC2),2)/length(fname)/12/4e-10))
title(['Spetrum-Average:', num2str(rotor_speed)])
saveas(h,['SpetrumAverage', num2str(rotor_speed),'.fig'])
save output_mat_CC1
toc


function DATA=V2Pa_Universal(Data,kulite_transform_ab)
    %first,check the size of Data and kulite_transform_ab,is same or not
    if size(Data,2)~=size(kulite_transform_ab,1)
        disp('转换参数不正确！！')
    end
     for k=1:size(Data,2)  
         DATA(:,k)= Data(:,k)*kulite_transform_ab(k,1)+kulite_transform_ab(k,2);%传感器B1
     end
       
end

 function [Pulse,rotor_speed]=keyRotation(thekeyData,fs)
   Tachometer = thekeyData;
    Threshold = 4;
    Temp_Num =  find(Tachometer>Threshold); 
    k = 0;
    for i = 1:1:length(Temp_Num)-1
        if (Temp_Num(i+1) - Temp_Num(i)) > 1
            k = k+1;
            the_Pulse(k,1) = Temp_Num(i);
        end
    end
    
    fr = fs/(mean(diff(the_Pulse)));
    the_rotor_speed = fr*60;
    Pulse=the_Pulse;
    rotor_speed=floor(the_rotor_speed/10)*10;    
 end
   