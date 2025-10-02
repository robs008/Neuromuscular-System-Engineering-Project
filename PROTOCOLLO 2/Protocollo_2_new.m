clear, close, clc;
warning off;

% I Plot sono stati tutti commentati per rendere più scorrevole il codice.

%% APERTURA DATI

addpath('sorgenti_gruppo7\')
data_open = readtable ('20_aperture.txt', 'Delimiter','\t');
data_close = readtable ('20_chiusure.txt', 'Delimiter','\t');
data_open_max = readtable ('max_apertura.txt', 'Delimiter','\t');
data_close_max = readtable ('max_chiusura.txt', 'Delimiter','\t');
data_test = readtable ('test_chiusura_apertura.txt', 'Delimiter','\t');

emg_close=table2array(data_close);   
emg_open=table2array(data_open);  
emg_close_max=table2array(data_close_max);   
emg_open_max=table2array(data_open_max);  
emg_test=table2array(data_test);  

% La prima colonna di questi dati rappresenta l'asse temporale,
% successivamente l'abbiamo tolto tanto la frequenza di campionamento è
% fissa a 2000Hz. Le altre 5 colonne rappresentano i 5 differenziali
% attorno alla parte più prossimale dell'avambraccio.

emg_test_new=emg_test(1:277205,2:6);

% % PLOT test
% figure
% for i=2:6
% 
%     plot((i-1)*4000+emg_test(:,i));
%     hold on;
% end
% 
% emg_test_new=emg_test(1:277266,2:6);
% 
% % PLOT test
% figure
% for i=1:5
% 
%     plot((i-1)*5000+emg_test_new(:,i));
%     hold on;
% end

label_test=[ones(22350,1);...
    zeros(37381-22350,1);...
    zeros(53640-37381,1);...
    zeros(69922-53640,1);...
    ones(87432-69922,1);...
    ones(102335-87432,1);...
    zeros(118718-102335,1);...
    ones(133363-118718,1);...
    ones(151513-133363,1);...
    ones(165118-151513,1);...
    zeros(181462-165118,1);...
    zeros(197889-181462,1);...
    ones(213951-197889,1);...
    zeros(230155-213951,1);...
    zeros(246348-230155,1);...
    ones(263304-246348,1);...
    ones(277205-263304,1)];
% figure
% plot(emg_test_new(:,2));
% hold on;
% plot(max(emg_test_new(:,2))*label_test,'r');
% hold off;


% % PLOT CHIUSURA
% figure
% for i=2:6
% 
%     plot((i-1)*1000+emg_close(:,i));
%     hold on;
% end


% Durante l'acquisizione ci siamo sgarrati 2 aperture quindi ne abbiamo
% aggiunto 2 muahahahahahaha;
% Di questo ne discuteremo tra di noi martedì
emg_open_new=emg_open(:,2:6);

% % PLOT APERTURA 
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_open_new(:,i));
%     hold on;
% end

% Durante l'acquisizione abbiamo preso 3 chiusure in più del previsto,
% quindi le abbiamo rimosse. A quando troppo a quando niente...
emg_close_new=emg_close(:,2:6);

% % PLOT CHIUSURA MAX
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_close_max(:,i));
%     hold on;
% end
% title('MAX close');
% 
% % PLOT APERTURA MAX
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_open_max(:,i));
%     hold on;
% end
% title('MAX open')

emg_close_max_new=emg_close_max(1719:11165,2:end);
emg_open_max_new=emg_open_max(661:9793,2:end);

% % PLOT CHIUSURA MAX
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_close_max_new(:,i));
%     hold on;
% end
% title('MAX close');
% 
% % PLOT APERTURA MAX
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_open_max_new(:,i));
%     hold on;
% end
% title('MAX open')
% pause()
% close all;

%% 1. Filtrare opportunamente i dati raccolti (se necessario). Si ricorda che i segnali EMG acquisiti sono di tipo singolo differenziale (SD). 

% Riportiamo esattamente quello che ci ha fatto fare il dottorando
% nell'esercitazione. Mette due filtri in cascata, perchè usare un
% passabanda direttamente non gli piaceva...
fsamp=2000;
fny=fsamp/2;
f2= 5; 
f1= 2; 
Ws= f1/fny; 
Rp=1; 
Rs= 20; 
Wp= f2/ fny;
[n,Ws] = cheb2ord(Wp,Ws, Rp, Rs);
[b,a] = cheby2(n, Rs, Ws, 'high');

% % visualizazione del filtro
% figure
% freqz(b,a,5000,fsamp);
% title('HPF to remove movement artifacts'); 
% è un filtro IIR poiché la fase non è lineare. 
% Quando la fase è lineare, i coefficienti 'a' e 'b' sono simmetrici o antisimmetrici


% filtro passa-basso a 350 Hz per togliere il rumore di alta frequenza (stesso procedimento di prima usando cheby2)
f2= 350; 
f1= f2 + 20;
Wp= f2/ fny;
Ws= f1/fny; Rp=1;
Rs= 20;[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);
[B,A] = cheby2(n, Rs, Ws, 'low'); 

% % visualizazione del filtro
% figure
% freqz(B,A,5000,fsamp)
% title('LPF to remove high frequency noise')


% applicazione dei filtri (confronto filter - filtfilt):
% non posso usare 'filter' in cascata se ho filtro IIR perché induce dei
% ritardi, perché può essere instabile. 
% Se ho un filtro FIR, non ha senso usare 'filter' in doppia passata,
% perché FIR ha fase lineare, la doppa passata la facciamo per eliminare la distorsione di fase (quindi nei filtri IIR).
emg_close_new=filtfilt(b,a, filtfilt(B,A, emg_close_new));
emg_open_new=filtfilt(b,a, filtfilt(B,A, emg_open_new));
emg_test_new=filtfilt(b,a, filtfilt(B,A, emg_test_new));

%% Esclusione riposo
[b,a]=butter(4,10/1000,"low");
% freqz(b,a,5000,fsamp);
emg_env_open=filtfilt(b,a,abs(emg_open_new));
emg_env_close=filtfilt(b,a,abs(emg_close_new));
emg_env_test=filtfilt(b,a,abs(emg_test_new));

% % apertura
% % Plot inviluppo
% figure
% for i=1:5
% 
%     plot((i)*1000+emg_env_open(:,i));
%     hold on;
% end
% figure
% plot(emg_env_open(:,2));
th_open=10;
ind_op=find(emg_env_open(:,2)>th_open);
% figure
% for i=1:5
% 
%     plot((i)*1000+emg_open_new(ind_op,i));
%     hold on;
% end
% 
% % chiusura
% % Plot inviluppo
% figure
% for i=1:5
% 
%     plot((i)*400+emg_env_close(:,i));
%     hold on;
% end
% figure
% plot(emg_env_close(:,3));
th_close=18;
ind_cl=find(emg_env_close(:,3)>th_close);
% figure
% for i=1:5
% 
%     plot((i)*1000+emg_close_new(ind_cl,i));
%     hold on;
% end

% test
% Plot inviluppo
% figure
% for i=1:5
% 
%     plot((i)*300+emg_env_test(:,i));
%     hold on;
% end
% figure
% plot(emg_env_test(:,3));
th_test=15;
ind_test=find(emg_env_test(:,3)>th_test);
% figure
% for i=1:5
% 
%     plot((i)*1000+emg_test_new(ind_test,i));
%     hold on;
% end

%% 2. Disporre i segnali ottenuti in una matrice per colonne. Il risultato atteso sarà una matrice m × n
% dove m indica i campioni e n i canali. Disporre prima le contrazioni associate all apertura della mano
% e poi quelle di chiusura. Fare in modo che la lunghezza temporale delle contrazioni incluse nella matrice
% sia la stessa.

% Per portare i dati alla stessa lunghezza abbiamo prima applicato
% resample, ma dava problemi, quindi abbiamo usato interp1 che con
% bioingegneria della riabilitazione ci ha salvato la vita <3
% Abbiamo resemplato a 240000 campioni perchè in teoria durante
% l'acquisizione ad ogni contrazione di 3s veniva seguita una pausa di
% altri 3 secondi; considerando una frequenza di campionamento di 2000Hz
% e avendo 20 contrazioni risultano un numero di campioni pari a: 
% (3)*20*2000=240000
nuovot = linspace(0,1,120000);
for i=1:5
    t = linspace(0,1, length(emg_open_new(ind_op,i)));
    emg_open_resampled(:,i)=interp1(t, emg_open_new(ind_op,i), nuovot,"spline");
    t = linspace(0,1, length(emg_close_new(ind_cl,i)));
    emg_close_resampled(:,i)=interp1(t, emg_close_new(ind_cl,i), nuovot,"spline");
end

% % (3)*17*2000=102000
% nuovot = linspace(0,1,102000);
% for i=1:5
%     t = linspace(0,1, length(ind_test));
%     emg_test_resampled(:,i)=interp1(t, emg_test_new(ind_test,i), nuovot,"spline");
% end

% % Plot dei dati ricampionati
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_close_resampled(:,i));
%     hold on;
% end
% title('close resempled'); hold off;
% 
% figure
% for i=1:5
% 
%     plot((i)*5000+emg_open_resampled(:,i));
%     hold on;
% end
% title('open resempled'); hold off;


%Concatenazione di aperture e chiusure
emg_op_cl=[emg_open_resampled;emg_close_resampled];



%% 3. Creare una copia della matrice al punto 2 e sostituire i segnali grezzi con gli inviluppi dei
% segnali EMG mediante raddrizzamento e filtraggio passabasso a 10 Hz, ottenendo
% nuovamente una matrice di dimensione m × n
[b,a]=butter(4,10/1000,"low");
%freqz(b,a,5000,fsamp);
emg_env=filtfilt(b,a,abs(emg_op_cl));
emg_test_new=emg_test_new(ind_test,:);
emg_env_test=filtfilt(b,a,abs(emg_test_new));
% % Plot inviluppo
% figure
% for i=1:5
% 
%     plot((i)*1000+emg_env(:,i));
%     hold on;
% end
% title('envelope'); hold off;

%% 4. Creare un vettore di label (0 e 1 oppure -1 e 1) di lunghezza m. Ogni label sarà associata alla
% tipologia di contrazione esaminata.
% classe 1= chiusura
% classe 0= apertura
label_class=[zeros(120000,1);ones(120000,1)];

train_data_env=emg_env;
train_data_grezzi=emg_op_cl;
label_train=label_class;
test_data_env=emg_env_test;
test_data_grezzi=emg_test_new;
label_test=label_test(ind_test,:);

%% 5. Applicare in Matlab un algoritmo di Support Vector Machine (SVM) per classificare i dati.
rng(123) 
% Get the number of samples 
num_samples = size(train_data_grezzi, 1); 
 
% Generate a random permutation of indices 
shuffled_indices = randperm(num_samples); 
 
% Shuffle the data and labels using the shuffled indices 
train_data_grezzi_s = train_data_grezzi(shuffled_indices, :); 
label_train_s = label_train(shuffled_indices, :); 
train_data_env_s = train_data_env(shuffled_indices, :); 
 
optimOptions = struct('AcquisitionFunctionName', 'expected-improvement-plus', ... 
                      'UseParallel', true, 'Repartition', true);  % Enable parallel computing if available 

Mdl_SVM_grezzi = fitcsvm(train_data_grezzi_s, label_train_s, ..., 
    'KernelFunction', 'rbf', ... 
    'Standardize',true, ... 
    'IterationLimit', 1e3, ... 
    'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ... 
    'HyperparameterOptimizationOptions', optimOptions); 
 
Mdl_SVM_ENV = fitcsvm(train_data_env_s, label_train_s, ..., 
    'KernelFunction', 'rbf', ... 
    'Standardize',true, ... 
    'IterationLimit', 1e3, ... 
    'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ... 
    'HyperparameterOptimizationOptions', optimOptions);
 
% VALUTAZIONE ACCURATEZZA SUL TRAIN 
% dati grezzi 
predicted_label_SVM_grezzi_TRAIN = predict(Mdl_SVM_grezzi, train_data_grezzi_s); 
accuracy_SVM_grezzi_TRAIN = sum(predicted_label_SVM_grezzi_TRAIN == label_train_s) / numel(label_train_s); 
% Inviluppi 
predicted_label_SVM_env_TRAIN = predict(Mdl_SVM_ENV, train_data_env_s); 
accuracy_SVM_env_TRAIN = sum(predicted_label_SVM_env_TRAIN == label_train_s) / numel(label_train_s); 
 
 
% VALUTAZIONE ACCURATEZZA 
% dati grezzi 
predicted_label_SVM_grezzi_TEST = predict(Mdl_SVM_grezzi, test_data_grezzi); 
accuracy_SVM_grezzi_TEST = sum(predicted_label_SVM_grezzi_TEST == label_test) / numel(label_test); 
% Inviluppi 
predicted_label_SVM_env_TEST = predict(Mdl_SVM_ENV, test_data_env); 
accuracy_SVM_env_TEST = sum(predicted_label_SVM_env_TEST == label_test) / numel(label_test);

%% 6. Applicare la Linear Discriminant Analysis (LDA) ai dati ottenuti ai punti 2 e 3.

% Creazione modello dati grezzi
Mdl_LDA_grezzi = fitcdiscr(train_data_grezzi, label_train,'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('UseParallel',true),...
    'HyperparameterOptimizationOptions',struct('Repartition',true));

% Creazione modello inviluppo
Mdl_LDA_env = fitcdiscr(train_data_env, label_train,'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('UseParallel',true),...
    'HyperparameterOptimizationOptions',struct('Repartition',true));


% VALUTAZIONE ACCURATEZZA SUL TRAIN
% dati grezzi
predicted_label_LDA_grezzi_TRAIN = predict(Mdl_LDA_grezzi, train_data_grezzi);
accuracy_LDA_grezzi_TRAIN = sum(predicted_label_LDA_grezzi_TRAIN == label_train) / numel(label_train);
% Inviluppi
predicted_label_LDA_env_TRAIN = predict(Mdl_LDA_env, train_data_env);
accuracy_LDA_env_TRAIN = sum(predicted_label_LDA_env_TRAIN == label_train) / numel(label_train);

% VALUTAZIONE ACCURATEZZA
% dati grezzi
predicted_label_LDA_grezzi = predict(Mdl_LDA_grezzi, test_data_grezzi);
accuracy_LDA_grezzi = sum(predicted_label_LDA_grezzi == label_test) / numel(label_test);
% Inviluppi
predicted_label_LDA_env = predict(Mdl_LDA_env, test_data_env);
accuracy_LDA_env = sum(predicted_label_LDA_env == label_test) / numel(label_test);

%% 7. Confrontare i risultati ottenuti usando i 2 metodi (SVM e LDA).

%% 8. Implementare un classificatore che sfrutta la misura di cosine-similarity per assegnare il
% movimento ad una specifica classe. Per farlo, è necessario considerare i dati delle contrazioni
% massimali per ottenere dei prototipi e quelli usati per creare la matrice al punto 2 come test.
% Partendo dalle contrazioni massimali durante apertura e chiusura della mano, stimare i valori
% di ARV (finestratura di 250 ms – overlap 50%) per ogni canale e salvarli. Mediare i valori di
% ARV così ottenuti per ciascun canale, ottenendo dei vettori di dimensione n che costituiscono
% i prototipi corrispondenti alle 2 contrazioni. Per identificare il tipo di contrazione da dati di
% test, valutare il vettore di ARV per tale contrazione e stimare la cosine-similarity con ciascun
% prototipo.

window_length=250e-3*fsamp;
overlap_int=window_length/2;
iteration=floor(length(emg_open_max_new)/overlap_int)-1;
for i=1:iteration
    start=(i-1)*window_length/2+1;
    stop=start+window_length-1;
    arv_gt_open(i,:)=mean(abs(emg_open_max_new(start:stop,:)));
    arv_gt_close(i,:)=mean(abs(emg_close_max_new(start:stop,:)));
end
arv_gt_open=mean(arv_gt_open);
arv_gt_close=mean(arv_gt_close);


iteration_test=floor(length(emg_op_cl)/overlap_int-1);
for i=1:iteration_test
    start=(i-1)*window_length/2+1;
    stop=start+window_length-1;
    arv_cycle=mean(abs(emg_op_cl(start:stop,:)));
    cos_op(i)=arv_gt_open*arv_cycle'/(norm(arv_gt_open)*norm(arv_cycle));
    cos_cl(i)=arv_gt_close*arv_cycle'/(norm(arv_gt_close)*norm(arv_cycle));
    if cos_op(i)>cos_cl(i)
        predicted_label_cos_train(i)=0;
    else
        predicted_label_cos_train(i)=1;
    end
end
label_cos_train=[zeros(1,479),ones(1,480)];
accuracy_COS_train = sum(predicted_label_cos_train == label_cos_train) / numel(label_cos_train);


iteration_test=floor(length(test_data_grezzi)/overlap_int-1);
for i=1:iteration_test
    start=(i-1)*window_length/2+1;
    stop=start+window_length-1;
    arv_cycle=mean(abs(test_data_grezzi(start:stop,:)));
    cos_op=arv_gt_open*arv_cycle'/(norm(arv_gt_open)*norm(arv_cycle));
    cos_cl=arv_gt_close*arv_cycle'/(norm(arv_gt_close)*norm(arv_cycle));
    if cos_op>cos_cl
        predicted_label_cos_test(i)=0;
    else
        predicted_label_cos_test(i)=1;
    end
end

label_cos_test=zeros(1,iteration_test);
for i=1:iteration_test
    start=(i-1)*window_length/2+1;
    stop=start+window_length-1;
    if sum(label_test(start:stop)) > window_length/2
        label_cos_test(i)=1;
    end
end
accuracy_COS_test = sum(predicted_label_cos_test == label_cos_test) / numel(label_cos_test);


