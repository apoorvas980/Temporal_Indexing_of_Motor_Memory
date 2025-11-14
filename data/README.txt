This folder contains the CSV files used in the analyses for the project Temporal Indexing of Motor Memory.
Each file includes trial-by-trial behavioral data such as hand-angle corrections, error directions, and timing information.
All files load directly into the MATLAB scripts found in the /Analysis directory.

Key Variables Used in the Analyses for computing Adaptation, Temporal Contextual Change and Generalization

handAng

The absolute hand angle (in degrees) that the participant produced on each trial wrt the line joining start circle and target.
This is the main behavioral measure used across experiments.

newAdaptationAng

The trial-to-trial change in hand angle. It captures how much participants corrected their movement on a given trial.
This relative measure is used for the adaptation and temporal-effects analyses.

rotDir

The error direction on the current trial, coded as:
–1 for clockwise (–15°) error
+1 for counterclockwise (+15°) error
Used for computing change in hand angle on the same trial.

newRotDir

The error direction from the previous trial, aligned, so that it matches the adaptation measured on the current trial.
It makes sure that when we look at trial-by-trial change in hand angle on trial t, it is paired with the error that actually caused it (trial t–1). This step does not modify the data. It is simply a clean bookkeeping alignment so the learning signal and its cause sit on the same row, making it easier to compute trial-by-trial adaptation

Other Variables

The datasets include additional columns.
These are retained for completeness but are not used directly in the main analyses.

Files in This Folder

twocontext_data.csv: dataset for the two-context task
threecontext_data.csv: dataset for the three-context task
control_data.csv: control dataset 


