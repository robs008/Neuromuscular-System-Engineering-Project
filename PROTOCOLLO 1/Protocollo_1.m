clear, close all, clc;
warning off;

%% CARICAMENTO DATI
addpath("dati_mat_7\");
addpath("Funzioni\");
load('bicipite_2.mat')
data.bicipite_2.mono=double(Data{1, 1});
load('bicipite_4.mat')
data.bicipite_4.mono=double(Data{1, 1});
load('bicipite_6.mat')
data.bicipite_6.mono=double(Data{1, 1});
load('bicipite_8.mat')
data.bicipite_8.mono=double(Data{1, 1});

load('tricipite_2.mat')
data.tricipite_2.mono=double(Data{1, 1});
load('tricipite_4.mat')
data.tricipite_4.mono=double(Data{1, 1});
load('tricipite_6.mat')
data.tricipite_6.mono=double(Data{1, 1});
load('tricipite_8.mat')
data.tricipite_8.mono=double(Data{1, 1});

% le colonne sono i canali, i primi 8 del tricipite e gli ultimi 8 bicipite
field_names = fieldnames(data);
fc=2048;
IED=5e-3;

%% 1. Filtrare opportunamente i dati raccolti dal bicipite (non considerando i primi 2 secondi,
% considerati come una fase di transitorio)
% noi facciamo anche già il tricipite

% calcolo dei campioni da escludere 
% escludo dall'analisi i primi 2 secondi
camp_esclusi=2*fc;
[cMA,cAR]=rico_mod(0.01,6,50,fc);
n=6;
wn=[10,350]/(fc/2);
[b,a]=butter(n,wn,'bandpass');

for i=1:numel(field_names)

    data.(field_names{i}).mono=filtfilt(cMA,cAR,filtfilt(b,a,data.(field_names{i}).mono(camp_esclusi:end,:)));
end

%% 2. Stimare i singoli e i doppi differenziali (SD e DD)

for i=1:numel(field_names)
    % SINGOLI DIFFERENZIALI
    % 14 colonne (canali) primi 7 tricipite, ultimi 7 bicipite
    for k=1:length(data.(field_names{i}).mono(1,:))/2-1
        data.(field_names{i}).sd(:,k)=data.(field_names{i}).mono(:,k+1)-data.(field_names{i}).mono(:,k);
    end
    for k=1:length(data.(field_names{i}).mono(1,:))/2-1
        data.(field_names{i}).sd(:,k+length(data.(field_names{i}).mono(1,:))/2-1)= ...
        data.(field_names{i}).mono(:,k+1+length(data.(field_names{i}).mono(1,:))/2)-...
        data.(field_names{i}).mono(:,k+length(data.(field_names{i}).mono(1,:))/2);
    end
    
    % DOPPI DIFFERENZIALI
    % 12 colonne (canali) primi 6 tricipite, ultimi 6 bicipite
    for k=1:length(data.(field_names{i}).mono(1,:))/2-2
        data.(field_names{i}).dd(:,k)=data.(field_names{i}).mono(:,k+2)-...
        2*data.(field_names{i}).mono(:,k+1)+data.(field_names{i}).mono(:,k);
    end
    for k=1:length(data.(field_names{i}).mono(1,:))/2-2
        data.(field_names{i}).dd(:,k+length(data.(field_names{i}).mono(1,:))/2-2)=...
        data.(field_names{i}).mono(:,k+2+length(data.(field_names{i}).mono(1,:))/2)-...
        2*data.(field_names{i}).mono(:,k+1+length(data.(field_names{i}).mono(1,:))/2)+...
        data.(field_names{i}).mono(:,k+length(data.(field_names{i}).mono(1,:))/2);
    end
end

%% 3. Stimare la densità spettrale di potenza (PSD) per le diverse contrazioni con il metodo più
% opportuno. Considerare epoche di 250 ms e le tre tipologie di segnali: monopolare, SD e DD.
% Quali modifiche subiscono gli spettri a seguito dell’applicazione del filtro spaziale? Mostrare
% i risultati tramite un subplot da 4 righe (contrazioni) e 3 colonne (tipologia di segnale).
epoch_length=round(250e-3*fc);
for i=1:numel(field_names)/2
    [P_mean.(field_names{i}).mono,f]=pwelch(data.(field_names{i}).mono(:,length(data.(field_names{i}).mono(1,:))/2+1:end)...
        -mean(data.(field_names{i}).mono(:,length(data.(field_names{i}).mono(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);
    P_mean.(field_names{i}).sd=pwelch(data.(field_names{i}).sd(:,length(data.(field_names{i}).sd(1,:))/2+1:end)...
        -mean(data.(field_names{i}).sd(:,length(data.(field_names{i}).sd(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);
    P_mean.(field_names{i}).dd=pwelch(data.(field_names{i}).dd(:,length(data.(field_names{i}).dd(1,:))/2+1:end)...
        -mean(data.(field_names{i}).dd(:,length(data.(field_names{i}).dd(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);
end
field_tipo= fieldnames(data.(field_names{i}));

% PLOT PSD

contrazioni={"2kg"; "4kg"; "6kg"; "8kg"};
Titolo_colonne={"Monopolari"; "Singoli differenziali"; "Doppi differenziali"};

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

for i = 1:numel(field_names)/2
    for k = 1:3
        subplot(numel(field_names)/2, 3, 3*(i-1) + k);
        plot(f, P_mean.(field_names{i}).(field_tipo{k}));
        
        % Aggiungi etichetta 'PSD (uV^2)' al primo grafico di ogni riga
        if k == 1
            ylabel('PSD (µV^2)');
        end
        
        % Rimuovi numeri dall'asse delle frequenze per tutte le righe tranne l'ultima
        if i < numel(field_names)/2
            set(gca, 'XTickLabel', []);
        end
        

        
        % Aggiungi la legenda
        if i==1
            legend_labels = strcat('ch_', num2str((1:numel(P_mean.(field_names{i}).(field_tipo{k}))).'));
            legend(legend_labels, 'Location', 'northeast');
        end
        % Aggiungi etichetta 'Frequenza (Hz)' solo nell'ultima riga
        if i == numel(field_names)/2
            xlabel('Frequenza (Hz)');
        end
        % Aggiungi titolo alla prima colonna
        if i == 1
            title(Titolo_colonne{k});
        end
                % Aggiungi testo in verticale e in grassetto accanto all'ultimo grafico di ogni riga
        if k == 3
            text(1.05, 0.5, contrazioni{i}, 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold');
        end
    end
end


%PLOT PER VALUTARE EVENTUALI ZONE DI INNERVAZIONE
for i=1:numel(field_names)/2
    figure
    title(field_names{i}); hold on;
    for k=8:14
        plot((k-1)*1000+data.(field_names{i}).sd(:,k));
        hold on
    end
end


%%
i = 3;
figure
for k = 1:7
    plot((k-1)*700 + data.(field_names{i}).sd(19300:19340, k+7));
    hold on
end
% Definizione delle stanghette di riferimento
x_start = 43; % Puoi modificare questo valore per spostare la posizione della stanghetta
y_start = 0; % Puoi modificare questo valore per spostare la posizione della stanghetta

x_length = 5; % Lunghezza della stanghetta orizzontale in campioni (0.1 secondi)
y_length = 500; % Lunghezza della stanghetta verticale in uV

% Aggiungi stanghetta orizzontale
line([x_start, x_start + x_length], [y_start, y_start], 'Color', 'k', 'LineWidth', 2);
text(x_start + x_length / 2, y_start -150, '2.5 ms', 'HorizontalAlignment', 'center');

% Aggiungi stanghetta verticale
line([x_start, x_start], [y_start, y_start + y_length], 'Color', 'k', 'LineWidth', 2);
text(x_start - 1, y_start , '500 uV', 'Rotation', 90, 'VerticalAlignment', 'middle');

hold off


% Aggiungi la legenda con i nomi dei canali
legend_labels = strcat('ch', num2str((1:7).'));
legend(legend_labels);

% Rimuovi i numeri dagli assi x e y
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

hold off


%% 4. Stimare la velocità di conduzione (CV) durante i diversi task mediante diverse tecniche
% descritte a lezione. Si consiglia di realizzare un diagramma a barre (con relativo errorbar) dove
% in ascissa sono ripotate informazioni circa l’intensità della contrazione. 

epoch_length=250e-3*fc;
n_epoch=floor(length(data.(field_names{i}).sd(:,1))/epoch_length);


for i=1:numel(field_names)/2
    n_epoch=floor(length(data.(field_names{i}).sd(:,1))/epoch_length);
    for k=1:n_epoch
        fatigue_parameters.cv.(field_names{i})(k)=mle_CV_est( data.(field_names{i}).sd(((k-1)*epoch_length+1):(k*epoch_length),...
            10:-1:8)',IED,fc); %solo i primi 2 canali
    end
end


for i=1:numel(field_names)/2
    x(i)=mean(fatigue_parameters.cv.(field_names{i})(1:10));
    error(i)=std(fatigue_parameters.cv.(field_names{i})(1:10));
end

% PLOT VELOCITA' DI CONDUZIONE
figure ()
b=bar(x','grouped');
hold on;
y = b(1).XEndPoints;
errorbar(y',x,error','k','linestyle','none');
hold off;
x_label_bar = {'2kg', '4kg', '6kg', '8kg'};
legend('Valori medi', 'Standard deviation')
xticklabels(x_label_bar);
ylabel('Velocità di conduzione (m/s)');

% Test statistico per verificare il trend

for i=1:3
    dati_condizione_1 = fatigue_parameters.cv.(field_names{i})(1:10);  
    dati_condizione_2 = fatigue_parameters.cv.(field_names{i+1})(1:10); 
    % Esegui il test Wilcoxon signed-rank
    [p_val_cv_mean(i), ~, ~] = signrank(dati_condizione_1, dati_condizione_2);
end

for i=1:4
    X = fatigue_parameters.cv.(field_names{i})(1:10);  
    % 2. Test di normalità
    % Test di Shapiro-Wilk
   % Calcolo della curtosi
    k = kurtosis(X);
    
    % Stampa del valore della curtosi
    fprintf('Curtosi dei dati: %.4f\n', k);
    
    % Verifica della normalità basata sulla curtosi
    if abs(k - 3) < 0.5 % Tollerenza di 0.5, può essere aggiustata in base ai tuoi dati
        kurt(i)=true;
    else
        kurt(i)=false;
    end
end



%% 5. Rappresentare quattro fatigue plot (uno per contrazione) dove saranno riportate informazioni
% sulla ampiezza del segnale (ARV e RMS), frequenza mediana (MDF), frequenza media
% (MNF) e velocità di conduzione (CV). Utilizzare epoche da 250 ms e i segnali SD.
for i=1:numel(field_names)/2
    n_epoch=floor(length(data.(field_names{i}).sd(:,1))/epoch_length);
    for k=1:n_epoch
         [fatigue_parameters.f_mean.(field_names{i})(k,:),fatigue_parameters.f_median.(field_names{i})(k,:),...
             fatigue_parameters.arv.(field_names{i})(k,:),fatigue_parameters.rms.(field_names{i})(k,:)]=...
             fatigue_plot(data.(field_names{i}).sd(((k-1)*epoch_length+1):(k*epoch_length),...
             length(data.(field_names{i}).sd(1,:))/2+1:length(data.(field_names{i}).sd(1,:))/2+1+2),fc);
    end
end
field_fatigue=fieldnames(fatigue_parameters);
for i=2:numel(field_fatigue)
    for k=1:numel(field_names)/2
        fatigue_parameters.mean.(field_fatigue{i}).(field_names{k})=...
            mean(fatigue_parameters.(field_fatigue{i}).(field_names{k})')';
    end
end
fatigue_parameters.mean.cv=fatigue_parameters.cv;

for i=1:numel(field_names)/2
    figure
    for k=2:numel(field_fatigue)
        subplot(3,2,k-1)
        plot(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}));
        xlabel('Epoca');
        title(strrep(field_fatigue{k}, '_', ' '));
    end 
    k=1;
    subplot(3,2,[5 6])
    plot(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}));
    xlabel('Epoca');
    title(strrep(field_fatigue{k}, '_', ' '));
    sgtitle([strrep(field_names{i}, '_', ' '), " kg"]);
end

%% 6. Realizzare cinque grafici (organizzati secondo subplot su una sola riga, con sovrapposti i
% risultati di ciascuna contrazione) dove per ciascuno sono riportate le rette interpolanti i dati
% per i descrittori EMG. In altre parole, sarà richiesto un grafico per l’ARV, uno per RMS, uno
% per MDF, un altro per MNF e infine uno per CV. Commentare i risultati ottenuti. 

y_ax={"Velocità (m/s)";"Frequenza (Hz)";"Frequenza (Hz)";"Ampiezza (µV)";"Ampiezza (µV)"};
for i=1:numel(field_names)/2
    figure
    for k=2:numel(field_fatigue)-1
        subplot(3,2,k-1)
        plot(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}),'.');
        hold on;
        p = polyfit(1:length(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i})),...
            fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}), 1);
        y=p(1)*(1:length(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i})))+p(2);
        plot(y);
        xlabel('Epoca');
        ylabel(y_ax{k});
        title(strrep(field_fatigue{k}, '_', ' '));
        hold off;
    end
    k=1;
    subplot(3,2,[5 6])
    plot(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}),'.');
    hold on;
    p = polyfit(1:length(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i})),...
        fatigue_parameters.mean.(field_fatigue{k}).(field_names{i}), 1);
    y=p(1)*(1:length(fatigue_parameters.mean.(field_fatigue{k}).(field_names{i})))+p(2);
    plot(y);
    xlabel('Epoca');
    ylabel(y_ax{k});
    title(strrep(field_fatigue{k}, '_', ' '));
    hold off;
    sgtitle(strrep(field_names{i}, '_', ' '));
end


%%

% 
% for i=1:numel(field_names)/2
%     figure
%     for k=2:numel(field_fatigue)
%         subplot(3,2,k-1)
%         data = fatigue_parameters.mean.(field_fatigue{k}).(field_names{i});
%         plot(data, '.');
%         hold on;
% 
%         % Calculate the linear fit
%         x = 1:length(data);
%         p = polyfit(x, data, 1);
%         y = polyval(p, x);
%         retta.(field_fatigue{k}).(field_names{i})=y;
% 
%         % Calculate R^2 value
%         ymean=mean(data);
%         SStot = sum((data - ymean).^2);
%         SSres = sum((y - ymean).^2);
%         rsq = SSres / SStot;
%         % Plot the linear fit
%         plot(x, y);
% 
%         % Add slope and R^2 text
%         text(mean(x), mean(data), ['Slope: ' num2str(p(1)) ', R^2: ' num2str(rsq)], 'FontSize', 8, 'BackgroundColor', 'w');
% 
%         xlabel('Epoca');
%         title(strrep(field_fatigue{k}, '_', ' '));
%         hold off;
%     end
% 
%     k = 1;
%     subplot(3,2,[5 6])
%     data = fatigue_parameters.mean.(field_fatigue{k}).(field_names{i});
%     plot(data, '.');
%     hold on;
% 
%     % Calculate the linear fit
%     x = 1:length(data);
%     p = polyfit(x, data, 1);
%     y = polyval(p, x);
%     retta.(field_fatigue{k}).(field_names{i})=y;
%     % Calculate R^2 value
%     ymean = mean(data);
%     SStot = sum((data - ymean).^2);
%     SSres = sum((y - ymean).^2);
%     rsq = SSres / SStot;
% 
%     % Plot the linear fit
%     plot(x, y);
% 
%     % Add slope and R^2 text
%     text(mean(x), mean(data), ['Slope: ' num2str(p(1)) ', R^2: ' num2str(rsq)], 'FontSize', 8, 'BackgroundColor', 'w');
% 
%     xlabel('Epoca');
%     title(strrep(field_fatigue{k}, '_', ' '));
%     hold off;
% 
%     sgtitle(strrep(field_names{i}, '_', ' '));
% end
% 
% colors = lines(numel(field_names)/2);  % Using MATLAB's default color order
% 
% figure
% for k=1:numel(field_fatigue)
%     subplot(5,1,k)
%     hold on;  % Ensure all plots are held on the same subplot
%     for i=1:numel(field_names)/2
%         % Normalize the retta values
%         retta.norm.(field_fatigue{k}).(field_names{i}) = retta.(field_fatigue{k}).(field_names{i}) / retta.(field_fatigue{k}).(field_names{i})(1);
% 
%         % Plot the normalized retta values with a specific color
%         plot(retta.norm.(field_fatigue{k}).(field_names{i}), 'Color', colors(i,:), 'DisplayName', strrep(field_names{i}, '_', ' '));
%         hold on;
%     end
%     hold off;  % Release the hold after all plots are done
%     xlabel('Epoca');
%     ylabel('Normalized Value');
%     title(strrep(field_fatigue{k}, '_', ' '));
%     if k==1
%         legend('show');  % Show legend with field_names
%     end
% end




%% 7. Filtrare opportunamente i dati acquisiti dal tricipite (non considerando i primi 2 secondi)
% Già fatto sopra

%% 8. Sommare i segnali monopolari del bicipite durante una contrazione e quelli del tricipite
% durante un'altra contrazione, in modo da simulare delle co-contrazioni.

% Useremo i dati dell'esercizio con 6kg.
% Ricampioniamo i dati del tricipite alla stessa lunghezza dei dati del
% bicibite.
% Nella somma considereremo il 100% del segnale del bicipite e il 60% del
% segnale del tricipite.

bic=data.bicipite_6.mono;
tri=data.bicipite_6.mono;
nuovot = linspace(0,1,length(data.bicipite_6.mono(:,1)));
for i=1:length(data.tricipite_6.mono(1,:))
    t = linspace(0,1, length(data.tricipite_6.mono(:,1)));
    tri(:,i)=interp1(t, data.tricipite_6.mono(:,i), nuovot,"spline");
end
mixed_data.mono=bic+tri;
% PLOT DEI SEGNALI
figure
subplot(2,2,1)
title('BICIPITE'); hold on;
for k=1:16
    plot((k-1)*3000+bic(:,k));
    hold on
end

subplot(2,2,2)
title('TRICIPITE'); hold on;
for k=1:16
    plot((k-1)*3000+tri(:,k));
    hold on
end

subplot(2,2,[3 4])
title('MIXED'); hold on;
for k=1:16
    plot((k-1)*3000+mixed_data.mono(:,k));
    hold on
end

%% 9. Ripetere i punti 2, 3, 4 in modo da valutare l’effetto del crosstalk (simulato) sui dati del
% muscolo target. 

% 9.2 STIMA SD E DD
% SINGOLI DIFFERENZIALI
% 14 colonne (canali) primi 7 tricipite, ultimi 7 bicipite
for k=1:7
    mixed_data.sd(:,k)=mixed_data.mono(:,k+1)-mixed_data.mono(:,k);
end
for k=9:15
    mixed_data.sd(:,k-1)=mixed_data.mono(:,k+1)-mixed_data.mono(:,k);
end

% DOPPI DIFFERENZIALI
% 12 colonne (canali) primi 6 tricipite, ultimi 6 bicipite
for k=1:6
    mixed_data.dd(:,k)=mixed_data.mono(:,k+2)-2*mixed_data.mono(:,k+1)+mixed_data.mono(:,k);
end
for k=9:14
    mixed_data.dd(:,k-2)=mixed_data.mono(:,k+2)-2*mixed_data.mono(:,k+1)+mixed_data.mono(:,k);
end

% PLOT DEI SEGNALI
figure
subplot(3,1,1)
title('Mixed mono'); hold on;
for k=1:16
    plot((k-1)*3000+mixed_data.mono(:,k));
    hold on
end

subplot(3,1,2)
title('Mixed sd'); hold on;
for k=1:14
    plot((k-1)*1000+mixed_data.sd(:,k));
    hold on
end

subplot(3,1,3)
title('Mixed dd'); hold on;
for k=1:12
    plot((k-1)*1000+mixed_data.dd(:,k));
    hold on
end

% 9.3 PSD
epoch_length=round(250e-3*fc);
[mixed_data.PSD.mono,f]=pwelch(mixed_data.mono(:,length(mixed_data.mono(1,:))/2+1:end)...
    -mean(mixed_data.mono(:,length(mixed_data.mono(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);
mixed_data.PSD.sd=pwelch(mixed_data.sd(:,length(mixed_data.sd(1,:))/2+1:end)...
    -mean(mixed_data.sd(:,length(mixed_data.sd(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);
mixed_data.PSD.dd=pwelch(mixed_data.dd(:,length(mixed_data.dd(1,:))/2+1:end)...
    -mean(mixed_data.dd(:,length(mixed_data.dd(1,:))/2+1:end)), bartlett(epoch_length), 0, fc, fc);

% PLOT PSD
figure
for k=1:3
    subplot(2,3,k)
    plot(f,P_mean.bicipite_6.(field_tipo{k}));
    title(strcat("No crosstalk ",field_tipo{k}))
    xlabel('Frequenza (Hz)')
    ylabel('PSD (µV^2)')
end
for k=1:3
    subplot(2,3,3+k)
    plot(f,mixed_data.PSD.(field_tipo{k}));
    title(strcat("Crosstalk ",field_tipo{k}))
    xlabel('Frequenza (Hz)')
    ylabel('PSD (µV^2)')
end

%%
% 9.4 CV
epoch_length=250e-3*fc;
n_epoch=floor(length(mixed_data.sd(:,1))/epoch_length);
for k=1:n_epoch
    mixed_data.cv(k)=abs(mle_CV_est( mixed_data.sd(((k-1)*epoch_length+1):(k*epoch_length),...
        8:10)',IED,fc)); %solo i primi 3 canali
end
for k=1:n_epoch
    unmixed_data.cv(k)=abs(mle_CV_est( data.bicipite_6.sd(((k-1)*epoch_length+1):(k*epoch_length),...
        8:10)',IED,fc)); %solo i primi 3 canali
end

%% verifica della coerenza tra le 2 tecniche

x = unmixed_data.cv;
y = mixed_data.cv;

mean_x= mean(x);
mean_y= mean(y);

% OLP regression
b = std(y)/std(x);
a = mean(y) - b*mean(x);
r = corrcoef(x,y); 
r = r(2);
B = finv(0.95,1,n_epoch-2)*(1-r^2)/(n_epoch-2);
ci.b = b * (sqrt(B+1) + [-sqrt(B) sqrt(B)]);
ci.a = mean(y) - ci.b([2 1])*mean(x);

figure("Color","w")
ax(1) = subplot(1,2,1);
plot(x,y,'Marker','o','MarkerSize',12,'MarkerEdgeColor','w','MarkerFaceColor','c','LineStyle','none')
line([min(x) max(x)],[min(x) max(x)]*b+a,'LineWidth',6)
line([min(x) max(x)],[min(x) max(x)],'Color','k')
title("OLP Regression")
ylabel("CV crosstalk; m/s")
xlabel("CV original; m/s")
text(ax(1),15,50,sprintf("CI slope: (%3.2f,%3.2f)\nCI int.: (%3.2f,%3.2f)",[ci.b ci.a]),"FontSize",12,"FontName","Arial")
legend('CV cross/CV original','Regressione lineare','Bisettrice');
% Bland-Altman plots
d.ls = x-y;
m.dls = mean(d.ls);
s.dls = std(d.ls);

ax(2) = subplot(1,2,2);
line((x+y)/2,x-y,'marker','o','linestyle','none','markersize',16,'MarkerEdgeColor','w','MarkerFaceColor','c')
line(ax(2).XLim,[1 1]*m.dls,'color',[1 1 1]*.7,'linewidth',4) % mean line
line(ax(2).XLim,[2 2]*s.dls+m.dls,'linestyle','--','color','k','linewidth',2) % +2SD line
line(ax(2).XLim,-[2 2]*s.dls+m.dls,'linestyle','--','color','k','linewidth',2) % -2SD line
title("Bland-Altman plot")
ylabel("\deltay (crosstalk - original; m/s)")
xlabel("mean CV (1/2(crosstalk + original ) m/s)")


set(ax,"FontName","Arial","FontSize",12)

% valutazione correlazione crosstalk vs no crosstalk
for i=1:4
    x=mixed_data.mono(:,8+i);
    y=bic(:,8+i);
    cos_sim.cross(i)=sum(x.*y)/(norm(x)*norm(y));
end

%% correzzione dello shift temporale
data_in_bss=[];
for k=1:n_epoch
    for i=1:4
        aus(:,i)=mixed_data.mono(((k-1)*epoch_length+1):(k*epoch_length),8+i);
        if i~=1
            shift=round(IED/mixed_data.cv(k)*fc);
            aus(:,i)=circshift(aus(:,i),i*shift);
        end
    end
    data_in_bss=[data_in_bss;aus];
end

%% PCA

[coeff_bicep, score_bicep, ~, ~, var_perc] = pca(data_in_bss);
threshold = 5;
% Identificare le componenti del bicipite da rimuovere
components_to_remove = (var_perc < threshold);
% Ricostruire i dati del bicipite senza le componenti correlate
score_bicep_cleaned = score_bicep;
score_bicep_cleaned(:, components_to_remove) = 0;
data_out_pca = score_bicep_cleaned * coeff_bicep';


% Ricostruzione schift temporale
data_out_pca_shift=[];
for k=1:n_epoch
    for i=1:4
        aus(:,i)=data_out_pca(((k-1)*epoch_length+1):(k*epoch_length),i);
        if i~=1
            shift=-round(IED/mixed_data.cv(k)*fc);
            aus(:,i)=circshift(aus(:,i),i*shift);
        end
    end
    data_out_pca_shift=[data_out_pca_shift;aus];
end



% valutazione tecnica PCA
for i=1:4
    x=mixed_data.mono(1:n_epoch*epoch_length,8+i);
    y=data_out_pca_shift(:,i);
    cos_sim.pca(i)=sum(x.*y)/(norm(x)*norm(y));
end

%% ICA
num_comp=4;
mdl_ica = rica(data_in_bss,num_comp);
score_bicep = transform (mdl_ica, data_in_bss);

for i=1:num_comp
    x=mixed_data.mono(1:n_epoch*epoch_length,1);
    y=score_bicep(:,i);
    cos_choose(i)=abs(sum(x.*y)/(norm(x)*norm(y)));
end

[~,ind]=max(cos_choose);

score_bicep(:,ind)=zeros(length(data_in_bss(:,1)),1);
% Reconstruction of the signals
data_out_ica=score_bicep*mdl_ica.TransformWeights';

% Ricostruzione schift temporale
data_out_ica_shift=[];
for k=1:n_epoch
    for i=1:4
        aus(:,i)=data_out_ica(((k-1)*epoch_length+1):(k*epoch_length),i);
        if i~=1
            shift=-round(IED/mixed_data.cv(k)*fc);
            aus(:,i)=circshift(aus(:,i),i*shift);
        end
    end
    data_out_ica_shift=[data_out_ica_shift;aus];
end



% valutazione tecnica PCA
for i=1:4
    x=mixed_data.mono(1:n_epoch*epoch_length,8+i);
    y=data_out_ica_shift(:,i);
    cos_sim.ica(i)=sum(x.*y)/(norm(x)*norm(y));
end




