
#Predicting RNA Secondary Structure by Genetic Algorithm.
##Modified Nussinov Algorithm

##Classes
* GeneticAlgorithm
* RNA
* Statistics
* Data

## Command Line
Basic Usage:
~~~~
python GeneticAlgorithm.py --help
~~~~

Run Test: test7 (there are four tests:test4,test5,test6,test7)
~~~~
python GeneticAlgorithm.py -t test7
~~~~

Drawing Structure (Must have Vienna Package Installed [3])
*Saves Postscript*
~~~~
python GeneticAlgorithm.py -t test4 -d
~~~~

Gathering Statistics Usage (Plot Mean/Max/Min Scores for each Generation)
~~~~
python Statistics.py --scatter
~~~~

## References
[1] Nussinov, Ruth, and Ann B. Jacobson. "Fast algorithm for predicting the secondary structure of single-stranded RNA." Proceedings of the National Academy of Sciences 77.11 (1980): 6309-6313.

[2] Belegundu, A. D., and T. R. Chandrupatla. "Optimization concepts and applications in engineering. 1999."

[3] Lorenz, Ronny, et al. "ViennaRNA Package 2.0." Algorithms for Molecular Biology 6.1 (2011): 1.

[4] Zuker, Michael. "Mfold web server for nucleic acid folding and hybridization prediction." Nucleic acids research 31.13 (2003): 3406-3415.

[5] Godfrey, Conrad. "RNA Secondary Structure Prediction: The
Co-transcriptional effect on RNA folding" Department of Statistics at Oxford, 2016. https://www.stats.ox.ac.uk/__data/assets/pdf_file/0008/6110/Godfrey.pdf




