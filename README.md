# Ingegneria del Sistema Neuromuscolare  
Progetto – A.A. 2023/2024  

Questo progetto si articola in due protocolli sperimentali basati sull’analisi di segnali elettromiografici (EMG).  

---

## 📌 Protocollo 1 – Classificatore mioelettrico per il controllo protesico  
- Obiettivo: Creare un classificatore basato su EMG dell’avambraccio per distinguere apertura e chiusura della mano.  
- Partecipanti: 1 soggetto sano.  
- Metodi: Acquisizione EMG da 5 muscoli dell’avambraccio (2 kHz), preprocessing, rettifica e inviluppo. Classificazione con SVM, LDA e cosine similarity.  
- Analisi: Confronto tra dati grezzi e inviluppi.  

### 🔬 Risultati principali
- LDA e SVM hanno raggiunto accuratezze medie intorno all’85% sugli inviluppi, superiori ai dati grezzi.  
- L’uso degli inviluppi ha migliorato la discriminazione riducendo il rumore e aumentando la stabilità del classificatore.  

---

## 📌 Protocollo 2 – Analisi della fatica muscolare e del crosstalk  
- Obiettivo: Studiare la fatica muscolare (bicipite e tricipite) e valutare strategie di riduzione del crosstalk.  
- Partecipanti: 1 soggetto sano + dati simulati.  
- Metodi: EMG ad alta densità (16 canali, 2048 Hz), estrazione di parametri di fatica (MDF, MNF, ARV, RMS, CV), confronto di configurazioni monopolari, singolo e doppio differenziale. Applicazione di PCA e ICA per la rimozione del crosstalk.  

### 🔬 Risultati principali
- Con l’aumento del carico si è osservata una diminuzione progressiva di MDF e MNF e un aumento di ARV e RMS, confermando l’andamento tipico dell’affaticamento muscolare.  
- La velocità di conduzione è aumentata con il carico (principio di Henneman).  
- I segnali differenziali (singolo e doppio) hanno ridotto il crosstalk in modo significativo rispetto ai segnali monopolari.  
- PCA ha rimosso il crosstalk con oltre il 97% di accuratezza, mentre ICA ha mostrato performance meno robuste.  

---

# Neuromuscular System Engineering  
Project – A.Y. 2023/2024  

This project is divided into two experimental protocols based on electromyographic (EMG) signal analysis.  

---

## 📌 Protocol 1 – Myoelectric classifier for prosthetic control  
- Aim: Develop a classifier based on forearm EMG to distinguish hand opening and closing.  
- Participants: 1 healthy subject.  
- Methods: EMG acquisition from 5 forearm muscles (2 kHz), preprocessing with Chebyshev filters, rectification and envelope extraction. Classification using SVM, LDA, and cosine similarity.  
- Analysis: Comparison between raw data and envelopes.  

### 🔬 Key Findings
- LDA and SVM achieved average accuracies around 85% on signal envelopes, outperforming raw data.  
- Cosine similarity yielded very high accuracy values (>95%), likely overestimated due to the metric.  
- Signal envelopes improved classification by reducing noise and increasing classifier stability.  

---

## 📌 Protocol 2 – Muscle fatigue and crosstalk analysis  
- Aim: Investigate muscle fatigue (biceps and triceps) and assess methods to reduce crosstalk.  
- Participants: 1 healthy subject + simulated data.  
- Methods: High-density EMG (16 channels, 2048 Hz), fatigue parameters (MDF, MNF, ARV, RMS, CV), comparison of monopolar, single differential, and double differential recordings. Application of PCA and ICA for crosstalk removal.  

### 🔬 Key Findings
- Increasing load led to a progressive decrease in MDF and MNF and an increase in ARV and RMS, consistent with typical muscle fatigue patterns.  
- Conduction velocity increased with load (Henneman’s principle), with anomalies at 8 kg suggesting technical or physiological limits.  
- Differential signals (single and double) significantly reduced crosstalk compared to monopolar recordings.  
- PCA removed crosstalk with over 97% accuracy, while ICA showed weaker robustness.  
