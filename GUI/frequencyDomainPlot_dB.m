function [the_freq,freq_dB]=frequencyDomainPlot_dB(signal,fs,Scale_data)
    
    data_fft = floor(fs/10);
    N_overlap = data_fft/2;
    N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
    data_freq = zeros(data_fft,1);
    for i = 1:N_seg
        data = signal(round((i-1)*(data_fft-N_overlap)+1):round((i-1)*(data_fft-N_overlap)+data_fft),:);
        temp_freq = abs(fft(data))*2/data_fft;
        data_freq = data_freq + temp_freq;
    end
    data_freq = data_freq/N_seg;
    freq_dB = 20*log10(data_freq/2e-5);     % 计算声压级
    the_freq = [1:round(data_fft/Scale_data) - 1]*fs/data_fft;  %数据频域离散刻度
    freq_dB=freq_dB(1:round(data_fft/Scale_data)- 1,:);
    plot(repmat(the_freq,1,size(signal,2)),freq_dB); 
 
%     set(gca,'YLim',[min(freq_dB)-1,max(freq_dB)+1])
%     ik = floor(fs/(Scale_data*rotor_speed(i_file)/60*29));
%     for nk=1:ik
%       plot(nk*rotor_speed(i_file)/60*29,min(freq_dB)-1,'^');
%       hold on
%     end
 end