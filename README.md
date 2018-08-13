# Mining NIH Research Abstracts and Determining Research Similarities
by Matt Stetz 2018

## Introduction
The sum_word_freq script mines the NIH database of published research abstracts and identifies keywords that define a researcher's corpus of published data. Keywords are identified by term frequency (TF).

Once the data are mined, you can run the cosine_calc script to determine the cosine similarity among all researcheres in a network

There are two scripts included in this repo:
* sum_word_freq3.py : mines NIH NCBI PubMed database
* cosine_calc.py : calculates cosine similarity for network analys

## Requirements
You need the following libraries:
* biopython
* elementtree
* numpy
* pylab/matplotlib
* seaborn

## Execution
The key input from the user is a txt file that lists the names of the researchers you wish to analyze in NCBI format. This should look like:

###### Smith J
###### Nguyen C
###### Lee JY
###### Jones MM
###### ...
###### Chen Z

You can also include additional NCBI search terms

###### Stetz MA AND University of Pennsylvania[ad]

There is also a list of trivial words that must be included 'trivial.txt' This file is similar to the english "stop words" used in sci-kit learn. I wrote all of the cleaning/munging functions myself for determining TFs so I included my own list of words to exclude.

## Example Results
Some example matrix plots are provided for various networks.
