# Melvin

<img src="https://github.com/XuemeiGu/MelvinPython/assets/37003667/7748ab20-3a26-4a3d-bd34-629b676d2b29" width="200"/>

[Automated Search for new Quantum Experiments](https://doi.org/10.1103/PhysRevLett.116.090405)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Phys. Rev. Lett. 116(9), 090405 (2016)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*Mario Krenn, Mehul Malik, Robert Fickler, Radek Lapkiewicz, and Anton Zeilinger*


Melvin is an innovative computer program designed to assist quantum physicists in the development and discovery of novel and beneficial quantum experiments. It was developed by [Mario Krenn](https://mariokrenn.wordpress.com/), marking a significant advancement in the field of quantum physics.

Here, I have partially reimplemented Melvin, originally developed in Mathematica, in Python utilizing SymPy for symbolic mathematics.

## Prerequisites:

To run the code from this repository, you need to install SymPy - **version 1.3**.

## Codes

Examples on creating quantum physical states and symbolic transformations: \
Related Mathematica codes are in the branch [Mathematica_Codes](https://github.com/XuemeiGu/MelvinPython/tree/Mathematica_Codes):

* The SimpleHOMExample program shows how to work with quantum states, and how the symbolic transformations work.
```
    SimpleHOMExample.nb
    SimpleHOMExample.py
```
* The CalcSRV program is a full version which searches for 3-particle high-dimensionally entanged states with existing optical elements.	
```
    CalcSRV.nb
    CalcSRV.py
```
* The SRVCaseCheck code is for checking whether a optical setup can produce 3-particle maximally entanged states. 
```
    SRVCaseCheck.py
```

* Here I add one more code for checking the SchmidtRankVector of a 4-particle maximally state\
  More details about Schmidt Rank Vectors, refer to [Structure of Multidimensional Entanglement in Multipartite Systems](https://doi.org/10.1103/PhysRevLett.110.030501). 
```
    Calc4ParticleSRV.nb
    Calc4ParticleSRV.py
```

## How to cite
if you want to cite the code for your work, you can use the following:
```
@misc{melvin_code,
  author = {Xuemei Gu},
  title = {Melvin: Automated Search for new Quantum Experiments -- Python Version},
  year = {2019},
  howpublished = {\url{https://github.com/XuemeiGu/MelvinPython}},
  note = {Accessed: 2019-04-03}
}
```

If you have any problems, please do not hesitate to contact me\
Email: njuxmgu@gmail.com or xuemei.gu@mpl.mpg.de
 
 

  
  


   
