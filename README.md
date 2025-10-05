# Automaticity-speeds-the-retrieval-of-instances-from-the-human-hippocampus
This dataset contains behavioral data, procedure code and analysis code with the manuscript "Automaticity speeds the retrieval of instances from the human hippocampus".
Data files (.mat) and associated scripts (.m) are divided into folders according to the subject of the analysis (e.g. Preprocess, TF, MVPA, Ripple etc.) .
All code is written in Matlab. 

General notes:
1)	Before running the analysis scripts, make sure you set the correct paths where the data file was extracted. All analyses were implemented in MATLAB.
2)	To run the code, the following open-source Matlab toolboxes are required:
•	Psychtoolbox (http://www.psychtoolbox.org), version: " Psychtoolbox 3". 
	D. H. Brainard, S. Vision, The psychophysics toolbox. Spat Vis 10, 433–436 (1997).
•	EEGLAB (https://sccn.ucsd.edu/eeglab/download.php), version: "eeglab2019". 
	A. Delorme, S. Makeig, EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. J. Neurosci. Methods. 134, 9-21 (2004).
3)  The dataset includes raw and analyzed data from Experiment 1, structured as follows:
•	BehaviorData: Contains behavioral data from the experiment.
•	RawData: Includes data used for SVM and ripple analyses, which contains folder ‘PowData’ and ‘RippleData’.
•	PowData: Contains theta band (4–8 Hz) activity recorded across electrode contacts in the hippocampus, middle-temporal gyrus, and prefrontal cortex for each participant.
•	RippleData: Includes raw iEEG signals bandpass-filtered in the 70–180 Hz range across hippocampus for each participant, other data of two brain areas (PFC and MTG) are not included due to huge size.
•	Fig1 – Fig4: Contain analysis results corresponding to each figure in the manuscript.
4) Note on Experiment 2
   Due to file size limitations, the data for Experiment 2 are not included in this repository. Please contact me directly if you require access to this portion of the data.

For further information, please contact: zhangyuanyuan_zoe@163.com


