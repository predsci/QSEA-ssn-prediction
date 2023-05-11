# QSEA-ssn-prediction
A novel quantile-based superposed-epoch analysis for predicting sunspot number

SSN prediction using Q-SEA model.

This code is intended to reproduce Figures 1-7 in the 
2023 solar physics paper: 
"On the Strength and Duration of Solar Cycle 25: 
A Novel Quantile-based Superposed-Epoch Analysis"
As such, it's not the most elegantly designed script. 
It's designed to generate each figure without having to 
change parameters interactively. 
An analysis version of this code provides more flexibility
but requires choosing/setting a number of parameters to 
get to the results in the paper. 

Please note that some of the later figures require that 
variables used in earlier figures be previously 
declared and populated with values. 

It can - and hopefully will - be further developed 
to be a more robust package. However, to be 
generally useful to the community, it should probably
be converted to a Python package. 

If you would like to use the logic/approach outlined 
in this code, please feel free to do so. I'd be 
grateful if you sent me an email (pete@predsci.com)
letting me know that you are using it, and please 
don't hesitate to contact me if you have any questions 
about the code. 

The code relies on two generally purpose libraries, 
both of which can be downloaded via CRAN. 
