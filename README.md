# adaptive-algorithm-for-irregular-Pareto-fronts
A Novel Archive Maintenance for Adapting Weight Vectors in Decomposition-   
based Multi-objective Evolutionary Algorithms, IEEE CEC 2021.

AMAWV.m is the algorithm. CDTLZ2 is an example of problem.  
Run the algorithm:     
main('-algorithm',Value,'-problem',Value,...) runs one algorithm on a problem with acceptable parameters.  
For example: main('-algorithm',@AMAWV,'-problem',CDTLZ2,'-N',200,'-M',3)  

All the acceptable parameters:  
%   '-N'            <positive integer>  population size  
%   '-M'            <positive integer>  number of objectives  
%   '-D'            <positive integer>  number of variables  
%	  '-algorithm'    <function handle>   algorithm function  
%	  '-problem'      <function handle>   problem function  
%	  '-evaluation'   <positive integer>  maximum number of evaluations  
%   '-run'          <positive integer>  run number  
%   '-save'         <integer>           number of saved populations  
%   '-outputFcn'	  <function handle>   function invoked after each generation  

