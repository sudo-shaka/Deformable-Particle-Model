To compile please make sure you have the required python depbendancies and C++ includes for pybind11 (Tested on Ubuntu). You may need to edit the path to header files if using another os. 

then:
python3 setupPyBind11.py build_ext --inplace

To use the model in python type:

import DPM

![out](https://user-images.githubusercontent.com/68864205/160499251-8eaa3176-4fd5-4e57-8e56-3761ebd4f069.gif)


Inspired by:
A. Boromand, A. Signoriello, F. Ye, C. S. O'Hern, and M. D. Shattuck, "Jamming of deformable polygons," Phys. Rev. Lett. 121 248003 (2018).


Requirements:
python3:
  - pybind11 : Get via "pip3 install pybind11"

C++ compilers:
  - build-essential : Get via "sudo apt install build-essential"
  - pybind11-dev : Get via "sudo apt install pybind11-dev"
