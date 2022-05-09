# README#

Ken Nakahara, 2022/03/10


This folder contains the programs about “Automatic Generation of Feedback Stabilizable State Space for Non-holonomic Mobile Robots” starting from April 2021 to March 2022.
All these programs are written by MATLAB (2021a).

Explanations about each folder are as follows.

* \practice

    This folder contains the codes for checking the control theory of the two-wheeled system.

* \sampling

    This folder contains the codes for collecting sensor variables and input samples for learning.

* \experiment_1d

    This folder contains the codes for mapping learning and feedback control when only the 1D coordinate transformation from the sensor variable s_3=θ to the virtual variable z2 (and the input transformations g and h) are unknown.
    Note that the adaptive grid distribution algorithm used in learning in this folder is unimproved compared with the algorithm used in the experiment_3d folder.

* \experiment_3d

    This folder contains the codes to for mapping learning and feedback control when the 3-dimensional coordinate transformations from the sensor variables s=(x,y,θ) to the virtual variables z=(z1,z2,z3) and the input transformations g and h are unknown.
    Note that the adaptive grid distribution algorithm used in the learning in this folder has been modified to achieve a three-dimensional coordinate transformations compared with the algorithm used in the experiment_1d folder.

* \experiment_3d_extend

    This folder contains the codes to for mapping learning and feedback control when the 3-dimensional coordinate transformations from the sensor variables s=(x,y,θ) to the virtual variables z=(z1,z2,z3) and the input transformations g and h are unknown.
    The codes included in this folder are intended to further extend the feedback controllable control area realized in the experiment_3d folder, and is currently under trial and error.
