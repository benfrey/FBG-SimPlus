# FBG-SimPlus V1.0

Fiber Bragg Grating (FBG) simulation tool for Finite Element Method (FEM) models. Features inclusion of temperature dependency and emulation within the program. The user can supply a data file and generate simulated reflection spectrums of an array of FBG sensors in response to:
* Longitudinal Strain (Uniform and Non-Uniform)
* Transverse Stress
* Temperature

![FBG-SimPlus Cover](python/GUI/resources/header.png)

### Quick Start Guide

To start FBG-SimPlus, execute the file "python/run.py". For further information on program usage, refer to the file "documentation.pdf" and the tutorial subdirectoy.

### Software Requirements

* Python Version: 3.8
* Required Modules: pyqt (Ver. 5), scipy, matplotlib, sympy, six, numpy
* Developed on macOS, supports Windows 10 and Linux

### Acknowledgements

B. Frey acknowledges support from the National Science Foundation and Department of Defense under Grants PHY-1659598 and PHY-1950744. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author and do not necessarily reflect the views of the National Science Foundation or Department of Defense.

### Considerations

Please email Ben Frey (freynben@gmail.com) with questions or with run-time issues.

Current known issues:
1. Instability when plotting spectral responses, sometimes leads to program crash.
2. Inability to exit from program on macOS, requires force quit.
3. Saving plot picture inconsistent functionality.

### Copyright

© 2020-2021 Benjamin N. Frey
