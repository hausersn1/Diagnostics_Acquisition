# Swept SFOAES

Stimulus frequency otoacoustic emissions are responses from the cochlea to pure tone stimuli. Since the stimulus is the same frequency as the emission, we have to use some tricks to separate the stimulus and the OAE. 

`Make_SFswept.m` generates a sweep of a probe tone and a suppressor, as well as sets some other stimulus parameters. This saves a `stim` struct so which is read by `Run_SFswept`

`Run_SFswept.m` generates the stimuli and plays them with the ER10X. The `stim` struct saves all the raw the data and the stimulus parameters. 

`Analyze_SFswept` will analyze a loaded stimuli structure that contains the data and plot the amplitude and phase. 

`rampsound.m` and `scalesound.m` are needed to make the stimuli. The ER10x files are needed for running this in the SNAP/ARDC lab. 

`simpleAnalysis.m` just does an fft of the response. 

## Discrete SFOAEs
 The code for discrete tone SFOAEs is in the discrete branch. Now that continuous sweep SFOAEs are working, we shouldn't really need this. 
