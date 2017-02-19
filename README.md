**spotgear**: Subset Profiling and Organizing Tools for Gel Electrophoresis Autoradiography in R
------
> An R Package for Fitting Bayesian Two-Dimensional Image Dewarping Models and 
Estimating Disease Subsets and Signatures

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

![](inst/example_figure/animation.gif)

**Figure 1**: Correction for smooth non-rigid gel deformation. Shown here are, for each gel batch, 19 serum lanes at 50 interior molecular weight grid points with large and small weights toward the left and right, respectively. Solid blue dots are detected peaks deviating from its true weight. Each detected peak is connected to a red triangle that represents the most likely actual molecular weight landmark. The black vertical curves together show the deformation, with each black vertical curve connecting locations with identical molecular weights. The curves are drawn for each grid point.

![](inst/example_figure/pwl_after_dewarping.pdf)

**Figure 2**: Visualization of the aligned gels. For each sample lane, given a set of peak-to-landmark alignment indicator ${Z_{gij}}_j$, the dewarping model maps the observed peaks to the landmark locations. Assuming at least one peak is detected, we then perform piecewise linear compression or stretching anchoring at matched landmarks along with two endpoint landmarks $\nu_0$ and $\nu_{L+1}$. The approximation is in general slightly different from the smooth function $\mathcal{S}$. One can view the Bayesian dewarping model based on landmarks as first estimating $S$ using downsampled landmarks to encourage nearby peaks to be aligned. Once the peak-to-landmark indicators are estimated, we use them for obtaining the piecewise linear approximation to $\mathcal{S}$. In this figure, we used the estimated the \textit{maximum a posteriori} estimates  $\widehat{Z_{gij}}$ to construct this approximation.



Maintainer:
--------------------------

Zhenke Wu (zhenkewu@umich.edu)

Department of Biostatistics

University of Michigan
