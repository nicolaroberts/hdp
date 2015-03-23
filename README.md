# hdp
R pkg for Hierarchical Dirichlet Process

To install:
```R
install.packages("devtools")
devtools::install_github("nicolaroberts/hdp")
```

R package to model categorical count data with a hierarchical tree of Dirichlet Processes. Includes functions to initialise a HDP of any shape, perform Gibbs sampling of the posterior distribution, and analyse the output. The underlying theory is described by Teh et al. (Hierarchical Dirichlet Processes, Journal of the American Statistical Association, 2006, 101:476). This R package was adapted from open source MATLAB and C code written by Yee Whye Teh and available here http://www.stats.ox.ac.uk/~teh/research/npbayes/npbayes-r21.tgz

```
Copyright (c) 2015 Genome Research Ltd. 
Author: Nicola Roberts <nr3@sanger.ac.uk> 
 
This program is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License version 3 
as published by the Free Software Foundation. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details <http://www.gnu.org/licenses/>. 
```

Copyright statement on original MATLAC and C code written by Yee Whye Teh, downloaded from 
http://www.stats.ox.ac.uk/~teh/research/npbayes/npbayes-r21.tgz

```
(C) Copyright 2004, Yee Whye Teh (ywteh -at- eecs -dot- berkeley -dot- edu)
http://www.cs.berkeley.edu/~ywteh

Permission is granted for anyone to copy, use, or modify these
programs and accompanying documents for purposes of research or
education, provided this copyright notice is retained, and note is
made of any changes that have been made.
 
These programs and documents are distributed without any warranty,
express or implied.  As the programs were written for research
purposes only, they have not been tested to the degree that would be
advisable in any important application.  All use of these programs is
```
