# Overview

Melvin is an algorithm to automatically design and find new useful experiments for quantum physicists, which was developed by [Dr. Mario Krenn](https://mariokrenn.wordpress.com/). The original codes are implemented in the software [Wolfram Mathematica](https://www.wolfram.com/mathematica/), of which the important features is that it can do symbolic, as well as numerical calculations.

Here I reimplement Melvin into [Python](https://www.python.org/) environment, which is open source software. [SymPy](https://www.sympy.org/en/index.html) is a Python library, which can be used for symbolic mathematics.


## Codes

Examples on creating quantum physical states and symbolic transformations:

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

## More infromation

For more information about Melvin, and other helpful links, take a look at these resources:

* **[Automated Search for new Quantum Experiments](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.090405)**,
  original paper for learning about Melvin
  
* **[Mathematcia codes](https://mariokrenn.wordpress.com/research/)**,
  in the mathematica_codes branch are two examples of Melvin provided by Dr. Mario Krenn.
  
 * **Experiments from Melvin’s solutions**, here I list some solutions from Melvin that have been realized in the laboratory of [Anton Zeilinger](https://www.iqoqi-vienna.at/people/zeilinger-group/anton-zeilinger/) at University Vienna & [IQOQI Vienna](https://www.iqoqi-vienna.at/research/zeilinger-group/quantum-entanglement-in-high-dimensional-systems/).
 
    1.**[Experimental Greenberger–Horne–Zeilinger entanglement beyond qubits](https://www.nature.com/articles/s41566-018-0257-6)**
    
    2.**[Multi-photon entanglement in high dimensions](https://www.nature.com/articles/nphoton.2016.12)**
    
    3.**[High-Dimensional Single-Photon Quantum Gates: Concepts and Experiments](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.180510)**
  
I appreciate all kinds of help, so thank you if you'd like to contribute. If you have any problems, please do not hesitate to contact me (email: njuxmgu@gmail.com).
 
 

  
  


   
