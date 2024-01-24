# IEMB-IPLF Filter
# Synopsis
This repository includes a MATLAB implementation of the Information exchange Multi-Bernoulli Filter [1] implemented using the Iterated Posterior Linearisation Filter (IPLF), the Extended Kalman Filter (EKF), the Iterated Kalman Filter (IEKF) and the Unscented Kalman Filter (UKF).
Also incuded is the Independant Multi-Bernoulli Filter which is implented using the IPLF and the UKF.
# Usage

Run scenario1.m to run all the above filters for Scenario 1 in the paper, which considers a Bernoulli birth model with large spatial uncertainty.

Run scenario2.m to run all the above filters for Scenario 2 in the paper, which considers a multi-Bernoulli birth model with each Bernoulli having a small spatial uncertainty.

If one is only interested in running one filter, the lines of code associated to the other filter can be commented.


# Contents

The function assign2D.m comes from [2]. 

The code of the generalised optimal subpattern assignment (GOSPA) metric has been obtained from [3].


[1] E. S. Davies and Á. F. García-Fernández, "Information Exchange Track-Before-Detect Multi-Bernoulli Filter for Superpositional Sensors," in IEEE Transactions on Signal Processing, vol. 72, pp. 607-621, 2024, doi: 10.1109/TSP.2024.3349769 

[2] D. F. Crouse, "The tracker component library: free routines for rapid prototyping," in IEEE Aerospace and Electronic Systems Magazine, vol. 32, no. 5, pp. 18-27, May 2017 https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary

[3] https://github.com/abusajana/GOSPA

