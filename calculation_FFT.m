%add comparison 2019-09-30 dgm
% i_method = 1, direct FFT without average; 
% i_method = 2, FFT with data segmentation and average; 

function [the_freq,data_freq]=calculation_FFT(signal,fs,i_method)
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    
    if i_method == 1
        data_fft = 2^floor(log2(length(signal)));
        data = signal(1:data_fft,:);
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %数据频域离散刻度
        data_freq = fft(data)*2/data_fft;
        data_freq = data_freq(1:data_fft/2.56,:); 
    end
    if i_method == 2
        data_fft = fs/10;
        N_overlap = data_fft/2;
        N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
        data_freq = zeros(data_fft,1);
        for i_seg = 1:N_seg
            data = signal((i_seg-1)*(data_fft-N_overlap)+1:(i_seg-1)*(data_fft-N_overlap)+data_fft,:);
            temp_freq = fft(data)*2/data_fft;
            data_freq = data_freq + temp_freq;
        end
        data_freq=data_freq(1:data_fft/2.56,:)/N_seg;
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %数据频域离散刻度
    end
    
 end