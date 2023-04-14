# anaesthesia-depth-control
EEG-based depth of anaesthesia control.


# Introduction
In addition to industry, the medical field has also seen a significant increase in the use of information technology combined with advanced process control methods. In this sense, much research is being done on computerised drug dosing in general anaesthesia, which belongs to the era of so-called personalised medicine.

The depth of anaesthesia is most commonly quantified by analysis of the electroencephalogram (EEG) using the bispectral index (BIS) and the patient state index (PSI). The goal of these measurements is to more accurately assess the depth of hypnosis, although other measurements, such as the flow of anaesthetic into the patient's body, can also be used and, in conjunction with appropriate pharmacokinetic and pharmacodynamic models, provide the anesthesiologist with additional information about the patient's condition.

The goal of computer-assisted drug delivery is to reduce the anesthesiologist's workload, increase safety, and shorten the patient's postoperative recovery time by reducing the dosage of the drug used.

The complete study (for the time being only in Slovenian) is also available as a pdf document.

Matlab scripts with Simulink models can be used to test different methods of controlling the depth of anaesthesia:

 - Open loop control
     - Target Controled Infusion (TCI)
     - Model Reference Control with Particle Swarm Optimisation (PSO)
- Closed loop control
    - PID regulator
    - Predictive Functional Controller (PFC)

# Associated files

## Requirements


[Open Source Patient Simulator](https://www.mathworks.com/matlabcentral/fileexchange/85208-open-source-patient-simulator) should be installed from MATLAB Central File Exchange and added to PATH.

| ![](/docs/img/simulator_concept.png) |
|:--:| 
| Patient Simulator concept – Clara Mihaela Ionescu (2023). [Open Source Patient Simulator](https://www.mathworks.com/matlabcentral/fileexchange/85208-open-source-patient-simulator), MATLAB Central File Exchange. |


## Matlab scripts

*The parameters of the controllers or models can be varied in scripts at predefined places to get different results*

- `main_pid.m` - [scipt](#pid) for testing PID controller.

- `main_tci.m` - [scipt](#tci) for testing TCI.

- `main_pso.m` - [scipt](#tci-pso) for testing TCI PSO.

- `main_pfc.m` - [scipt](#pfc) for testing PFC.

- `pfc_vs_tci.m` - [scipt](#pfc-vs-tci) for automated side-by-side PFC vs TCI comparison.

- `main_pfc_pd_analysis.m` - [scipt](#pfc-pd-analysis) for automated PFC analysis with wrong estimated PD model parameters.

- `main_pfc_detail_analysis.m` - [scipt](#pfc-pd-detail-analysis) for automated PFC analysis with wrong estimated PK and PD model parameters.

## Simulink models

- `Patient_PID_noDeadTime.slx` - model used for PID control without taking dead time into account.

- `Patient_PID.slx` - model used for PID control.

- `Patient.slx` - model used for PFC (main topic).

- `Patient_TCI.slx` - model used for TCI.

- `Patient_TCI_PSO.slx` - model used for TCI PSO.

## Additional files
- `Q_rand_<X>pct.mat` - matrix used for simulating wrong estimated parameters of PK model (X% range of system matrix)


# Demo

## PID control<a id='pid'></a>

User can use model that takes (long) dead time into account. If that dead time is not taken into account, control algorithm is not stable.
```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Select PID type  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim('Patient_PID_noDeadTime');
%sim('Patient_PID');
```
PID control **without** taking dead time into account:
| ![](/docs/img/reg_bis_pid_noDeadTime.png) |
|:--:| 
| PID without dead time. |

PID control **with** taking dead time into account:
| ![](/docs/img/reg_bis_pid.png) |
|:--:| 
| PID with dead time. |

## Target Controled Infusion<a id='tci'></a>

TCI with **exact** estimated parameters:
| ![](/docs/img/tci_exactParam_noDist.png) |
|:--:| 
| TCI with exact parameters. |


TCI with **wrong** estimated parameters of PK model (10% range of system matrix). We can see bias:
| ![](/docs/img/tci_10pctParam_noDist.png) |
|:--:| 
| TCI with wrong estimated parameters. |

## Target Controled Infusion with PSO<a id='tci-pso'></a>

TCI PSO with **exact** estimated parameters:
| ![](/docs/img/tci_pso_exactParam_noDist.png) |
|:--:| 
| TCI PSO with exact parameters. |


## Predictive Functional Controller<a id='pfc'></a>

PFC with **exact** estimated parameters:
| ![](/docs/img/reg_bis_pfc.png) |
|:--:| 
| PFC with exact parameters. |


## Side-by-side PFC vs TCI comparison<a id='pfc-vs-tci'></a>

PFC vs TCI comparison with **exact** estimated parameters:
| ![](/docs/img/PFC_vs_TCI_exactParam_noDist.png) |
|:--:| 
| PFC vs TCI with exact parameters. |

PFC vs TCI comparison with with **wrong** estimated parameters of PK model (20% range of system matrix):
| ![](/docs/img/PFC_vs_TCI_20pctParam_noDist.png) |
|:--:| 
| PFC vs TCI with wrong estimated parameters. |

## PFC analysis<a id='pfc-pd-analysis'></a>

Sensitivity of parameter $\sigma$ analysis:
| ![](/docs/img/PFC_pd_analysis_sigma.png) |
|:--:| 
| Parameter $\sigma$ analysis. |

Sensitivity of parameter $\gamma$ analysis:
| ![](/docs/img/PFC_pd_analysis_gamma.png) |
|:--:| 
| Parameter $\gamma$ analysis. |

Sensitivity of parameter $T_d$ (dead time estimation) analysis:
| ![](/docs/img/PFC_pd_analysis_mrtviCas.png) |
|:--:| 
| Parameter $T_d$ analysis. |

Sensitivity of parameter $C_{50P}$ analysis:
| ![](/docs/img/PFC_pd_analysis_c50p.png) |
|:--:| 
| Parameter $C_{50P}$ analysis. |

Sensitivity of parameter $C_{50R}$ analysis:
| ![](/docs/img/PFC_pd_analysis_c50r.png) |
|:--:| 
| Parameter $C_{50P}$ analysis. |

## PFC detailed analysis<a id='pfc-pd-detail-analysis'></a>

In the simulations used to obtain the results presented above to analyse the impact of parameter misestimation, we have assumed that we know the exact parameters of the PK models. However, if we do not know them precisely, the final results may of course differ. In the case that we misestimate the parameters of the PK model in the range of 20\% and we also misestimate the parameters of the PD model, the simulation results can be seen:
| ![](/docs/img/PFC_pk_pd_analysis.png) |
|:--:| 
| PFC detailed analysis. |


