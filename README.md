# README_program #

Ken Nakahara, 2022/03/10


## English ##

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


## 日本語 ##

このフォルダには，2021年4月から2022年3月まで行った卒業研究「非ホロノミック移動ロボットのフィードバック安定化可能な状態空間の自律的生成」のプログラムがあります．
これらのプログラムは全てMATLAB (2021a) で書かれています．

各フォルダについての説明は以下の通りです．

* \practice

    このフォルダには，二輪車両の制御理論の確認のためのコードが入っています．

* \sampling

    このフォルダには，学習のためのセンサ変数および入力のサンプル収集に関するコードが入っています．

* \experiment_1d

    このフォルダには，センサ変数s_3=θから仮想変数z2への1次元の座標変換（および入力変換g,h）のみを未知とした場合における写像学習，およびフィードバック制御を行うためのコードが入っています．
    なお，このフォルダ内の学習で使われている適応的格子分布アルゴリズムは，experiment_3dフォルダで使われているアルゴリズムと比べて未改良であることに注意してください．

* \experiment_3d

    このフォルダには，センサ変数s=(x,y,θ)から仮想変数z=(z1,z2,z3)への3次元の座標変換および入力変換g,hを未知とした場合における写像学習，およびフィードバック制御を行うためのコードが入っています．
    なお，このフォルダ内の学習で使われている適応的格子分布アルゴリズムは，experiment_1dフォルダで使われているアルゴリズムと比べて3次元の座標変換を実現するために改良されています．

* \experiment_3d_extend

    このフォルダには，センサ変数s=(x,y,θ)から仮想変数z=(z1,z2,z3)への3次元の座標変換，および入力変換g,hを未知とした場合における写像学習，およびフィードバック制御を行うためのコードが入っています．
    このフォルダが含むコードは，experiment_3dフォルダで実現されているフィードバック制御可能な制御範囲を更に広げることを目的としており，現在試行錯誤しているコードです．






