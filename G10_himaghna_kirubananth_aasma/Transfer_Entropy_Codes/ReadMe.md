###   Transfer Entropy Implementation

All the Helper Functions and The main implementation of the Transfer Entropy (TE)
is given in the files.

The Script's results that were shown in the Presentation were the results of *main_discrete.m*

(Could not complete the with other methods of TE estimation as the code was not well optimized)

(found some online libraries and toolboxes for the calculation of TE but I was unable to resolve compatibility issues)
(The online resources were: TRENTOOL - a MATLAB toolbox for TE estimation with many implementation methods

                            IDTxl - Information Dynamics Library for Python with JIDT implementation (JAVA based) 
                            (Explore IDTxl as python implementation in MATLAB was possible by using *pyenv* function
                            -sad that I was not able to make use of it early on)
                            )

In the Script if you want to use other methods you have to change the TE estimation function from *discrete_pipepline* to 
either *kraskov_pipeline* or *conditional_pipeline*

Along with the code the the matrices of TE and the figures of the results are also provided with.
