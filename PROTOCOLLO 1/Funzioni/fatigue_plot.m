function [f_mean_p,f_median_p,arv_p,rms_p]=fatigue_plot(x,fsamp)
%
% Funzione per calcolare i parametri relativi ai fatigue plot
%
% Input:    x:          Segnale
%
% Output:   f_mean_p:   Frequenza media
%           f_median_p: Frequenza media
%           arv_p:      Average rectified value
%           rms_p:      Root mean square
%
% Autori:   Roberto Pilotto
%           Salvatore Rapisarda
%
%
%
%
% Periodogramma semplice
epoch_len=length(x(:,1));
[P,f] = pwelch(x-mean(x), bartlett(epoch_len), 0, epoch_len, fsamp); 

% Calcolo frequenza mediana
A=0.5*sum(P);
for i=1:length(A)
    Ak=0;
    k=1;
    while Ak<A
        Ak=Ak+P(k,i);
        k=k+1;
    end
    f_median_p(i)=f(k);

    % Calcolo frequenza media
    f_mean_p(i) = sum(P(:,i).*f)/ sum(P(:,i));
    
    % Calcolo ARV
    arv_p(i)=mean(abs(x(:,i)));
    
    % Calcolo RMS
    rms_p(i)=rms(x(:,i));
end



end


