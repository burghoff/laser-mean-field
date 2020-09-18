Simulate a QCL (or any other laser) using the mean-field theory discussed in Burghoff's "Frequency-modulated combs as phase solitons" (ArXiv:2006.12397, 2020). This theory reduces the Maxwell-Bloch equations of a laser down to a single equation that can be integrated over the round trip of a cavity, by normalizing the internal field to the internal gain of the laser above threshold. 

Specifically, it uses a symmetric split-step method to find the F function in equation (5). Also produces a theoretical plot for the theoretical form of the soliton (equation (9), currently only valid for a Fabry-Perot cavity with either R1=1 or R2=1).

Input: a parameter structure generated using NLFM_Params

Output: a solution structure that stores the field evolution

Example call:

`s=Mean_Field_Theory(NLFM_Params('kpp',-1000))`

(Sets all parameters to their default but the dispersion, which is -1000 fs^2/mm.)

Contributors: David Burghoff, Levi Humbard
