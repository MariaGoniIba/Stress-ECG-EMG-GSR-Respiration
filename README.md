# Feature selection in mental stress analysis using multiple biological signals
The aim of this project is to classify different levels of mental stress using biological signals such as electrocardiogram (ECG), heart rate (HR), galvanic skin response (GSR), electromyogram (EMG) and respiration.

This is the code for the paper "Feature selection in mental stress analysis using multiple biological signals".

## ABOUT DATASET
### Online Repository Link
[Link to data repository](https://physionet.org/content/drivedb/1.0.0/)

### Description of the data
This database, contributed to PhysioNet by its creator, Jennifer Healey, contains a collection of multiparameter recordings from healthy volunteers, taken while they were driving on a prescribed route including city streets and highways in and around Boston, Massachusetts. The objective of the study for which these data were collected was to investigate the feasibility of automated recognition of stress on the basis of the recorded signals, which include ECG, EMG (right trapezius), GSR (galvanic skin resistance) measured on the hand and foot, and respiration.

## METHODOLOGY
### Methods and aims
The aforementioned database contains a collection of biological signals corresponding to 18 volunteers. Those volunteers drive through a route along Boston’s streets. The route alternates periods of rest, highway and city, simulating three levels of stress. We consider the rest state as the lower level of stress, driving on highway as the medium level, and the periods along the city as the highest level of stress.

This project consists on a classifier of 3 different levels of stress. We first extract a large number of feature for each biological signal according to an extensive literature review. The feature selection method and classifier consists of a genetic algorithm and a least squares linear discriminant.

### Code
* *MainScript.m*: this script runs the entire project with all biological signals and prints the most selected features. It is possible to run this code selecting subsets of signals and/or changing the parameters of the genetic algorithm.
* *FeaturesECG.m*: function that creates ECG features extracting the RR interval with the function *intervalRR.m*.
* *FeaturesHR.m*: function that extract features from heart rate recordings
* *FeaturesGSR.m*: funtion that extract features from galvanic skin responses. Two modalities: hand or foot. It calls the function *fGSR.m* to extract impulse magnitudes, durations and number of responses.
* *FeaturesEMG.m*: function that extract features from the EMG signals. 
* *FeaturesRESP.m*: function that extract features from respiration signals. It calls the function *respsignal.m* that calculates the respiratory interval, detecting inhaling and exhaling times.
* *Classification.m*: it runs the classification procedure, consisting in a genetic algorithm calling the function *GeneticAlgorithm.m* and a Least Squares Linear Discriminant. To increase robustness of the model, we run the whole procedure n_rep = 20 times.

The flowchart below shows the followed procedure along with the scripts names.

<p align="center">
  <img src="https://github.com/MariaGoniIba/Stress-ECG-EMG-GSR-Respiration/blob/main/flowchart.png"
</p>

## RESULTS
  
Results for different features, time slots and parameters of the genetic algorithm can be found in the paper "Feature selection in mental stress analysis using multiple biological signals".
  
## CITATION
If you use this code, please consider citing our work:

Goñi, M., Mohino, I., Llerena, C., Gil-Pita, R. and Rosa, M. (2013) Feature selection in mental stress analysis using multiple biological signals. 10th IASTED International Conference on Signal Processing, Pattern Recognition and Applications. Innsbruck (Austria)

## PAPERS
* [Feature selection in mental stress analysis using multiple biological signals](https://www.actapress.com/Abstract.aspx?paperId=455013)
* [Detecting stress during real-world driving tasks using physiological sensors](https://ieeexplore.ieee.org/document/1438384)
* [Wearable and automotive systems for affect recognition from physiology](http://dspace.mit.edu/handle/1721.1/9067)
* [PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals](https://www.ahajournals.org/doi/full/10.1161/01.CIR.101.23.e215)
