# MOTS2
Multi-Objective Tabu Search 2

# file structure

For each application, the folders are:

- bin: The output executables go here, both for the app and for any tests and experiments.
- build: This folder contains all object files, and is removed on a clean.
- config: create the development and production config files
- doc: Any notes, like  assembly notes and configuration files, are here.
- include: All project header files. All necessary third-party header files that do not exist under /usr/local/include are also placed here.
- lib: Any libs that get compiled by the project, third party or any needed in development. Prior to deployment, third party libraries get moved to /usr/local/lib where they belong, leaving the project clean enough to compile on our Linux deployment servers. This is used to test different library versions than the standard.
- experiments: smaller classes or files to test technologies or ideas, and keep them around for future reference. They go here, where they do not dilute the real application’s files, but can still be found later.
- src: The application and only the application’s source files.
- test: All test code files. 

# python extension

python 3

# DevOps

Coverage:
[![Coverage Status](https://coveralls.io/repos/github/christos-tsotskas/MOTS2/badge.svg?branch=master)](https://coveralls.io/github/christos-tsotskas/MOTS2?branch=master)

# LICENCE

By cloning/downloading/forking this repository (or part of it), you automatically agree and accept the license terms.

Reference one of the following:

```
@article{tsotskas2015multi,
  title={Multi-Objective Tabu Search 2: First Technical Report},
  author={Tsotskas, Christos and Kipouros, Timoleon and Savill, Mark A},
  year={2015},
  publisher={Cranfield University}
}

@inproceedings{tsotskas2013biobjective,
  title={Biobjective optimisation of preliminary aircraft trajectories},
  author={Tsotskas, Christos and Kipouros, Timoleon and Savill, Mark},
  booktitle={International Conference on Evolutionary Multi-Criterion Optimization},
  pages={741--755},
  year={2013},
  organization={Springer Berlin Heidelberg}
}

@article{razzaq2013multi,
  title={Multi-objective optimization of a fluid structure interaction benchmarking},
  author={Razzaq, Mudassar and Tsotskas, C and Turek, S and Kipouros, T and Savill, M and Hron, J},
  journal={CMES: Computer Modeling in Engineering \& Sciences},
  volume={90},
  number={4},
  pages={303--337},
  year={2013}
}

@article{razzaq2010insight,
  title={Insight into Fluid Structure Interaction Benchmarking through Multi-Objective Optimization},
  author={Razzaq, M and Tsotskas, C and Turek, S and Hron, J and Kipouros, T and Savill, M},
  journal={Wiley InterScience},
  volume={2},
  pages={1--25},
  year={2010}
}


```

## license dependencies

MOTS2, on its own is licensed under the Apache model.

In addition, the following packages were also used that follow individual licensing schemes.

- easylogging++ , v9.96.4, https://github.com/muflihun/easyloggingpp, MIT license
