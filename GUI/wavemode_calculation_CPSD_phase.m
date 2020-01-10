%i_method = 1: 时间谱
%i_method = 2: 阶比谱
function [GAMMA,freq]=wavemode_calculation_CPSD_phase(data,fs,nk,phase,rotor_speed,save_directory,i_file,fname)
    %[Pulse,rotor_speed]=keyRotation(data(:,end),fs);
    %filename_label = ['CPSD method ', num2str(i_method)];
    signal = data(:,1:nk);
        
    L_signal = length(signal);
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
    for k=1:nk
        for l = 1:nk
            [C{k,l},freq] = cpsd(signal(:,k),signal(:,l),Wind,Noverlap,Nfft,fs);          
        end
    end
    GAMMA = zeros(Nfft/2+1,nk+1);
    mode=-nk/2:nk/2;
    for m = 1:nk+1
        temp_f = zeros(Nfft/2+1,1);
        for k = 1:nk
            for l = 1:nk%2*pi*k/nk
                temp_f = temp_f + 0.5*C{k,l}*exp(i*mode(m)*phase(k))*exp(-i*mode(m)*phase(l));
            end
        end
        GAMMA(:,m) = temp_f/(nk*nk);
    end
    
    GAMMA = 10*log10(abs(GAMMA)/4e-10);

    colormap(jet);
    imagesc([-nk/2:nk/2],freq,GAMMA); colorbar;
     axis xy; %ylim([1,30]);
     testTime='试验-2019-12-13';
    title({[testTime,'-截止模态分析',' -CPSD method '];[char(fname(i_file)),'-转速: ',num2str(round(rotor_speed))]},'FontSize',14)

    xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);colorbar;
%     

 end
   
