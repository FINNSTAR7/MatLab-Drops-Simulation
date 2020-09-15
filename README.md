# MatLab-Drops-Simulation

This project was started to verify the expected and nearly guaranteed number of runs a Warframe player needs to collect all of a Warframe's components at least one time each.


simul.m will count the number of trials it takes to succeed on each probability inputted for each sample. For example, inputting "`simul`" into the command line will promt the user to input the number of samples to use, the probabilites, and the parts. If the user inputs "`1000000`", "`0.5 0.5`", and "`1 1`" respectively, simul.m will count the number of trials until each sample succeeds at getting at least 1 of each part, where each part has a 50% chance to drop.

It will then calculate the mean, median, mode, 99th percentile, 99.9th percentile, and the 99.99th percentile from this data, and plot these values along with the data and multiple confidence intervals (99%, 95%, 90%, 75%, and 50%).


On Command Line:<br />
`simul`

On Promt (defaults):<br />
Number of Samples: `500000` (default)<br />
Probabilities: `''` (default)<br />
Parts: `1 1 1 ... 1` (default, array of ones the same length as the number of probabilties given)

On Promt:<br />
Number of Samples: `1000000`<br />
Probabilities: `1/6 1/3 1/2` (notice elements are spearated by spaces, not commas)<br />
Parts: `1 1 2` (this specifies that you want at least one of the first two parts, and at least two of the last part)

The promt will remember your last inputs as well, so you do not need to input everything again if you run simul again (in the current MatLab session).


Note that increasing the samples will increase accuracy, however it will also increase run time. The opposite is also true, decreasing the samples will decrease accuracy, but also decrease run time. Keep this in mind, as some computers may need a smaller sample size than others.


For any questions or additional assistance please contact me through the wiki here: https://warframe.fandom.com/wiki/Message_Wall:FINNER
