function plot_results_demo(filename)
    load(filename);
    time_eeg = (1:1:length(raw_eeg))/fs-1/fs;
    
    figure, set(gcf,'color','white'); 
    hold on;
    plot(time_eeg, raw_eeg(142,:),'b','LineWidth',3);
    hold on;
    plot(time_eeg, bcg_artifact(142,:),'r-','LineWidth',2);
    ylabel('EEG Amplitude');
    xlabel('Time (Seconds)');
    title('Dynamic Harmonic Regression (EGI Evaluation Data)');
    legend('Raw EEG', 'BCG Artifact'); legend('boxoff');
    
    
    figure, set(gcf,'color','white'); 
    hold on;
    plot(time_eeg, raw_eeg(142,:),'b','LineWidth',3);
    hold on;
    plot(time_eeg, clean_eeg(142,:),'g-','LineWidth',2);
    ylabel('EEG Amplitude');
    xlabel('Time (Seconds)');
    title('Dynamic Harmonic Regression (EGI Evaluation Data)');
    legend('Raw EEG', 'Clean EEG'); legend('boxoff');
end