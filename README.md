# MatLab-Drops-Simulation

This project was started to verify the expected and nearly guaranteed number of runs a Warframe player needs to collect all of a Warframe's components at least one time each.


simul.m will count the number of trials it takes to succeed on each probability inputted for each sample. For example, inputting "`simul(1/3, 1/3 ,1/3);`" into the command line will have simul.m count the number of trials until an individual sample succeeds on the first 1/3, the second 1/3, and third 1/3 probability at least once each. Once all probabilities have succeeded it will move on to the next sample, until the data for all 1,000,000 samples have been collected.

It will then calculate the mean, median, mode, 99th percentile, 99.9th percentile, and the 99.99th percentile from this data, and plot these values along with the data and multiple confidence intervals (99%, 95%, 90%, 75%, and 50%).


Examples of Command Line Usage:<br />
`simul(0.5, 0.5);` (coin flip)<br />
`simul(1/6, 1/6, 1/6, 1/6, 1/6, 1/6);` (d6 toss)<br />
`simul(.75);`<br />
`simul(10, 20, 30, 40);`


To alter the number of iterations simul() runs through (default is 1,000,000) edit the "`samples`" line directly under the "`Begin Simulation & Plot`" header in the simul() function itself.

Note that increasing the samples will increase accuracy, however it will also increase run time by approximetely the same factor samples was increased by. The opposite is also true, decreasing the samples will decrease accuracy, but also decrease run time by approximetely the same factor the samples were decreased by. Keep this in mind when altering the function for yourself, as some computers may need a smaller sample size than others.

For example, 10,000 samples may take 20 seconds. If increased by a factor of 50, for 500,000 samples, the run time will increase to approximetely 1,000 seconds.


For any questions or additional assistance please contact me through the wiki here: https://warframe.fandom.com/wiki/Message_Wall:FINNER
