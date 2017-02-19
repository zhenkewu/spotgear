**spotgear**: Subset Profiling and Organizing Tools for Gel Electrophoresis Autoradiography in R
------
> An R Package for Fitting Bayesian Two-Dimensional Image Dewarping Models

Details
-------------------------------------

Autoimmune diseases, e.g., scleroderma, are human immune system’s
    responses to autoantigens in which the body produces specific autoantibodies
    that target these autoantigens but also cause tissue damage. The autoantibody
    composition is strikingly different among patients, as indicated by the many
    different radiolabeled patterns obtained from the mass-spectrometry-based
    technology - gel electrophoresis autoradiograms (GEA). Human recognition of
    patterns is not optimal when the patterns are composite/complex, or patterns
    are scattered across a large number of sam- ples. However, multiple sources of
    error (including irrelevant intensity differences across gels and warping of
    the gels) have traditionally precluded automation of pattern discovery using
    autoradiograms. In this package, we overcome these limitations by novel initial
    gel data preprocessing (Bayesian two dimensional image dewarping/registration) and then provide methods to
    detect disease subgroups.

Animation of Fitted Bayesian Dewarping Result
-------------------------------

![**Figure**: Correction for smooth non-rigid gel deformation. Shown here are, for each gel batch, $19$ serum lanes at $50$ interior molecular weight grid points with large and small weights toward the left and right, respectively. Solid <span style="color:blue">blue</span> dots "<span style="color:blue">$\bullet$</span>" are detected peaks deviating from its true weight. Each detected peak is connected to a <span style="color:red">red</span> triangle "<span style="color:red">$\Delta$</span>" that represents the most likely actual molecular weight landmark. The black vertical curves together show the deformation, with each black vertical curve connecting locations with identical molecular weights. The curves are drawn for each grid point.](inst/example_figure/animation.gif)


Maintainer:
--------------------------

Zhenke Wu (zhenkewu@umich.edu)

Department of Biostatistics

University of Michigan
