# Exploring multi-functional desirability

In this Shiny app you can change the importance given to the 24 ecosystem functions used in the synthesis analysis. The easiest way to run this App is via opening an R console and typing:

```R
library(shiny)

runGitHub("treeweb_synthesis","TerrestrialEcology-ugent",subdir="shiny")
```

This might take a while before the Shiny app appears since the rather large prediction data files (~50MB) will be downloaded into a temporary directory.

In the App you can: (i) look at the default results, (ii) apply equal importance to all functions or (iii) manually set the importance. The graphs can be saved on the computer via the save button. Do not hesitate to mail me graphs that you find interesting or if you have further questions.

If you expect to run the App multiple time and you are bored of having to wait for the downloads of the data files each separate time you can do the following:

* clone the repository into your local computer (via git clone https://github.com/lionel68/multi_desire.git)

* in the R console set the working directory to the folder with the cloned files and do: runApp()
