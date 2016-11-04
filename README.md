Data analysis of the \Geyser\ data set
=========================================

This repository contains the handed in version of an end of term
project in computational statistics. I took this course during the
winter semester 2015/16 at Heidelberg University. The aim was to use
tools developed in the lecture to analyse a simple data set. I decided
to use the Geyser data set, which ships with every R installation.

I tried to build a framework in order to find good predictive models
for the 'Old Faithful Geyser'. Therefore I introduced a skew cost
function to account for the asymmetry in the problem (it is worse, if
someone misses an eruption than if the has to wait a little bit). With
this I had to work out some fit routines as the standard regression
tools do onlye suppert the standard, quadratic cost function. In the
end I tried to measure the predictive capabilities of my prediction
candidates.
