
A federated learning algorithm to quantify racial disparities in kidney graft failure rates using US registry data from 29,468 patients across 149 transplant centers
==============================================


## Outline
1. Description
2. dGEM-disparity workflow
3. Data analysis results


## Description
In the United States, there exist racial disparities in kidney transplant access and post-transplant outcomes between Black and White patients. Site of care is considered a key contributor to the disparities. Using multi-site data to examine the effect of site of care on racial disparities, hospital-specific effects are essential and usually estimated with a generalized linear mixed effect model. However, due to regulations for protecting patients’ privacy, patient-level data typically cannot be shared across centers/hospitals. As a result, we developed a federated learning framework, dGEM-disparity (decentralized algorithm for Generalized linear mixed Effect Model for disparity quantification). Consisting of two parts, dGEM-disparity first provides accurately estimated common effects and calibrated hospital-specific effects by requiring only aggregated data from each center, and then adopts a counterfactual modeling approach to assess whether the graft failure rates differ if Black patients had been admitted at transplant centers in the same distribution as White patients were admitted. With the Organ Procurement and Transplantation Network (OPTN) registry data collected on 35,497 patients from 200 transplant centers, we found that if Black patients had been admitted to the centers following the distribution of White patients, 29 fewer per 10,000 Black patients (95% CI: [27, 32]) would die or have graft failure within one year after receiving a kidney transplant on average. Our proposed dGEM-disparity is generally applicable to other disparities research related to differential access to care. 


## dGEM-disparity workflow 
<img width="700" alt="image" src="https://user-images.githubusercontent.com/38872447/181802378-8a1499ce-0490-421e-8fc5-d6563da25543.png">




Workflow of the proposed decentralized algorithm for generalized mixed effect models framework for disparity quantification -- dGEM-disparity. This framework includes two main parts. Part I: the dGEM algorithm includes the initialization of the common effects (i.e., association between outcome and covariates) (Step I) and the estimation of center-effects (Step II) with the meta-analysis estimates of the common effects, following with a center-level calibration via shrinkage estimator. Part II, the counterfactual rates for disparity quantification are calculated with the estimates from Part I. The details of Step III in Part II are presented in the following section.


## Data analysis results
<img width="700" alt="image" src="https://user-images.githubusercontent.com/38872447/181802626-548d5cc2-b295-4c11-95a5-8b493564fa88.png">


The observed graft failure rate for the Black patients was 8.60%. The estimated counterfactual graft failure rate decreased by 0.29% (95% CI, [0.27%, 0.32%]), down to 8.31%. This means that if Black patients had been admitted to the centers following the distribution of White patients, 29 fewer per 10,000 Black patients (95% CI: [27, 32]) would have died or had graft failure within one year after receiving a kidney transplant on average.
