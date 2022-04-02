To compile please make sure you have the required python depbendancies and C++ includes for pybind11 (Tested on Ubuntu). You may need to edit the path to header files if using another os. 

then:
python3 setupPyBind11.py build_ext --inplace

To use the model in python type:

import DPM

![out](https://user-images.githubusercontent.com/68864205/161360831-56d26cb1-8585-4685-a922-a52f455e2364.gif)


Inspired by:
A. Boromand, A. Signoriello, F. Ye, C. S. O'Hern, and M. D. Shattuck, "Jamming of deformable polygons," Phys. Rev. Lett. 121 248003 (2018).


Requirements:
python3:
  - pybind11 : Get via "pip3 install pybind11"
  - setuptools : Get via "pip3 install setuptools" or "sudo apt install python3-setuptools"

C++ compilers:
  - build-essential : Get via "sudo apt install build-essential"
  - pybind11-dev : Get via "sudo apt install pybind11-dev"
