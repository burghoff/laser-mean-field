Simulate a QCL (or any other laser) using the mean-field theory discussed in Burghoff's "Unraveling the origin of frequency modulated
combs using active cavity mean-field theory" (https://doi.org/10.1364/OPTICA.408917, preprint version [here](http://arxiv.org/abs/2006.12397)). This theory reduces the Maxwell-Bloch equations of a laser down to a single equation that can be integrated over the round trip of a cavity, which is then iterated. A GUI allows certain parameters to be adjusted on the fly if desired.

Specifically, it uses a symmetric split-step method to evolve the F function according to equation (5). It also produces a theoretical plot for the theoretical form of the soliton (equation (9), currently only valid for a Fabry-Perot cavity with either R1=1 or R2=1).

Input: a parameter structure generated using NLFM_Params

Output: a solution structure that stores the field evolution

| Call | Explanation |
|---|---|
|`s=Mean_Field_Theory(NLFM_Params('kpp',-1000))`|Sets all parameters to their default but the dispersion, which is -1000 fs^2/mm.|
|`s=Mean_Field_Theory(NLFM_Params('gc',0,'useLinPwrApprox',1))`|Disables gain curvature and makes the linear power approximation: the result quickly converges to the phase version of the NLSE.|
|`s=Mean_Field_Theory(NLFM_Params('initphi','cosine2','numTr',5000))`|Initialize the phase to a periodic cosine, which converges to the N=2 harmonic state. (For numerical reasons it eventually decays to the fundamental.)|
|`s=Master_Equations(NLFM_Params('numTr',5000,'ampphase',1));`|Solves the master equations in amplitude-phase mode, comparing the results to the extendon theory.|

Contributors: David Burghoff, Levi Humbard
