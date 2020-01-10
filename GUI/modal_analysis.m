function [nk,GAMMA] = modal_analysis(DATA,Fs,circle1,save_directory,fname3)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
%    DATA=Data(:,mic_transform_ab(:,2));
   Phase=circle1/180*pi;
  DATA=  DATA./rms(  DATA);%对信号做正则化处理,正则化稀疏度更好

   [the_freq1,freq_dB1]=frequencyDomainPlot_dB2(  DATA,Fs,2.56);
    nk=32;
    
    BPF=the_freq1(find(freq_dB1(:,1)==max(freq_dB1(:,1))));rotor_speed=BPF/17;
    listName={'1BPF';'2BPF';};
    [GAMMA,freq]=wavemode_calculation_CPSD_phase(  DATA,Fs,nk,Phase,rotor_speed*60,save_directory,1,fname3);
    df=freq(2)-freq(1);list=[round(16.5*rotor_speed/df) round(17.5*rotor_speed/df);round(33.5*rotor_speed/df) round(34.5*rotor_speed/df);];

    for k=1:2
    [x(k),y(k)]=find(GAMMA==max(max(GAMMA(list(k,1):list(k,2),:))),1);Zdata=GAMMA(x(k),:);Zdata(y(k))=0;[y_1(k)]=find(Zdata==max(Zdata),1);
    o_RI(k)=x(k)*df;m_RI(k)=y(k)-(nk/2+1);m_RI_1(k)=y_1(k)-(nk/2+1);
    text(m_RI(k),o_RI(k),{[num2str(round(m_RI(k)*10)/10),',',num2str(round(GAMMA(x(k),y(k))))];[num2str(o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60))]});text(m_RI_1(k),o_RI(k),[num2str(m_RI_1(k)),',',num2str(round(GAMMA(x(k),y_1(k))))]);   
    end
    for k=1:2
    
    bar([-nk/2:nk/2],GAMMA(x(k),:));hold on;
    ylim([40 70]);
    title([char(fname3(i_file)),'-',num2str(round(rotor_speed*60)),'rpm','-',listName{k}],'FontSize',14)
    saveas(h1{k},[save_directory,'\',strrep(char(fname3{k}),'.wav','-'),'modeplot-',listName{k},'.png'])
    saveas(h1{k},[save_directory,'\',strrep(char(fname3{k}),'.wav','-'),'modeplot-',listName{k},'.fig'])
    end
end

