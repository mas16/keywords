# Scrape, Clean, and Mine NIH Research Abstracts and Determine Research Similarities
by MAS 2018

## Introduction
The ```wordfreq2.py``` script scrapes the NIH database of published research abstracts and identifies keywords that define a researcher's corpus of published data. Keywords are identified by comparing term frequency (TF) values.

Once TFs are known, you can run the ```cosine_calc.py``` script to determine the cosine similarity among all researcheres in a network

There are two scripts included in this repo:
* ```wordfreq2.py``` : scrapes and mines NIH NCBI PubMed database
* ```cosinecalc.py``` : calculates cosine similarity for network analysis

## Requirements
You need the following libraries:
* biopython
* elementtree
* numpy
* pylab/matplotlib
* seaborn

## Execution
The key input from the user is a ```.txt``` file that lists the names of the researchers you wish to analyze in NCBI format. This should look like:

> ###### <p> Smith J </p>
> ###### <p> Nguyen C </p>
> ###### <p> Lee JY </p>
> ###### <p> Jones MM </p>
> ###### <p> ... </p>
> ###### <p> Chen Z </p>

You can also include additional NCBI search terms, for example:

> ###### <p> Stetz MA AND University of Pennsylvania[ad] </p>

There is also a list of trivial words that must be included, ```trivial.txt``` . This file is similar to the english "stop words" used in scikit-learn. I wrote all of the cleaning/munging functions myself for determining TFs so I included my own list of words to exclude.

First fill out the the user-specified data directories in wordfreq then run that script. After that run ```cosine_calc.py``` (you do not need to specifiy anything).

## Example Results
Two files are generated by wordfreq for each author analyzed:

* A bar plot showing the author's highest TFs (default is 10 highest TFs). An example is shown below:

![](./JHU_TEST/BARRICK_D.png)

* An example correlation matrix for an entire network of scientists is shown below (network is the Biophysics Department at Johns Hopkins):

![](./JHU_TEST/matrix.png)

* The entire output from the Johns Hopkins example is provided in the ```JHU_TEST``` directory.
