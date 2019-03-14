# CAM2FIB
Interpretation and interpolation of G-code for creating stream files for FIB/SEM instruments.
The currently implemented output format is designed for use with an FEI Helios NanoLab G3 instrument, but may easily be adapted to other instruments.

This code accepts most G-code files as an input and was tested with the G-code dialects RS274D, WinCNC and FreeCAD.
In its current implementation the overall patterning time is equally distributed over all patterning points.
Machining or printing of 3D-structures is therefore currently limited to layered structures.

Please read the manual for further instructions on installing and running the interpreter program.
