# Overview

Melvin is an algorithm to design and find new useful experiments for quantum physicists, which was developed by [Dr. Mario Krenn](https://mariokrenn.wordpress.com/). 

The original codes are implemented in the software [Wolfram Mathematica](https://www.wolfram.com/mathematica/), which supports symbolic language very well. 

Here I reimplement Melvin into Python environment, which is open source software. [SymPy](https://www.sympy.org/en/index.html) is a Python library, which can be used for symbolic mathematics.
 
Links to Sections:

* [Installation](#installation)
* [Creating your first Quantum Program](#creating-your-first-quantum-program)
* [More Information](#more-information)


## Codes

Examples on creating quantum physical states and symbolic transformations:
```
    SimpleHOMExample.nb
    SimpleHOMExample.py
```
These program show how to work with quantum states, and how the symbolic transformations work.
	
```
    CalcSRV.nb
    CalcSRV.py
```
* The Mathematica code is a full version which searches for 3-particle high-dimensionally entanged states with existing optical elements.
* The python version shows how to check whether a optical setup can produce 3-particle maximally entanged states. (Currently, the example doesn't contain the searching experiments part.)
 


## More infromation

For more information about Melvin, and other helpful links, take a look at these resources:

* **[Automated Search for new Quantum Experiments](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.090405)**,
  original paper for learning about Melvin
  
* **[mathematcia codes](https://mariokrenn.wordpress.com/research/)**,
  two mathematcia examples of Melvin (SimpleHOMExample.nb, CalcSRV.nb), which are provided by Dr.Mario Krenn.

* **[another python codes](https://github.com/StephenCzy/Melvin_python_version)**,
  another python version of Melvin, which does not include criterias for checking maximun SRV and other complicated situations.   

  
  


   
