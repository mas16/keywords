"""
Calculate cosine similarity for bag of words.
Words are ranked by raw TF.
TFs from abstracts scraped from NCBI using sum_word_freq.py.

Note: the vector space is re-defined for each pairwise comparison.

By MAS 05/2018.

Updated by MAS 10/2018 - style clean up.
"""

###
#Import Libraries
import wordfreq
import matplotlib
#Need to change matplotlib backend for plotting
matplotlib.use('Agg')
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

#Title for the similarity matrix plot
TITLE = 'Similarity Scores for JHU PMB Faculty\n'

###
#Define constants

#Working directory path
BASEPATH = wordfreq.OUTPATH

#List of PIs path
LISTPATH = wordfreq.INPATH

#Path to output matrix plot
OUTPATH = wordfreq.OUTPATH+'matrix.png'

###
#Define functions

def get_listlen(listpath):
    """Get list of authors from text."""
    with open(listpath) as authors:
        data = authors.read()
        data = data.split('\n')[:]
    return data

def get_authors(basepath, authors):
    """Get author keywords from text."""
    authorpath = basepath+authors+'.txt'
    keywords = get_keywords(authorpath)
    return keywords

def organize_keywords(keywords):
    """Split keyword input from txt."""
    final = [entry.split() for entry in keywords]
    return final

def convert_float(v):
    """Convert TFs to floats to enable cosine calculation.

    Recall: data is read from txt so all objects are str.
    """
    v_list = []
    for x,entry in enumerate(v):
        if x % 2 == 0:
            v_list.append(entry)
        else:
            v_list.append(float(entry))
    return v_list
                   
def clean_keywords(keywords):
    """Treat author/search term line separately.

    This is necessary because parsing of first line -
    will vary depending on number of search terms.
    Author term is always first so just keep that one.
    """
    final = []
    for x, entry in enumerate(keywords):
        if x == 0:
            final.append([entry[0]])
        else:
            final.append(convert_float(entry))
    return final

def get_keywords(datapath):
    """Execute all keyword cleaing functions."""
    with open(datapath) as keywords:
        wordlist = keywords.readlines()
        wordlist = organize_keywords(wordlist)[:]
        wordlist = clean_keywords(wordlist)[:]
    return wordlist       

def construct_vectorspace(v1, v2):
    """Construct vector space for one pair of authors."""
    t1 = []
    for entry1 in v1:
        if entry1[0] not in t1:
            t1.append(entry1[0])
    for entry2 in v2:
        if entry2[0] not in t1:
            t1.append(entry2[0])
    return t1

def populate_vectorspace(v1, combined):
    """Populate text vector space with TFs from each author."""
    t1 = []
    for entry in combined:
        #Initialize counter
        c = 0
        for item in v1:
            if entry == item[0]:
                t1.append(item)
                c += 1
        #For specific case where TF is 0
        if c == 0:
            t1.append([entry,0])
    return t1

def square(v):
    """Square each element in list.

    Could probably just use an array here...
    """
    v2 = [entry[1]**2 for entry in v]
    return v2

def calc_cos(v1, v2):
    """Calculate cosine similarity."""
    prod = [v1[x][1]*v2[x][1] for x in range(len(v1))]
    v1_sq = square(v1)[:]
    v2_sq = square(v2)[:]
    cos = sum(prod)/(np.sqrt(sum(v1_sq))*np.sqrt(sum(v2_sq)))
    return cos

def plot_matrix(title, names, data, outpath):
    """Plot using seaborn heatmap."""
    ax = sns.heatmap(data, cmap='coolwarm',
                     xticklabels=names, yticklabels=names)
    ax.tick_params(direction='out', labelsize=8)
    plt.title(title, fontweight='bold')
    plt.tight_layout()
    plt.savefig(outpath)
    print ('\nCorrelation Matrix Saved')
    return 0

def cos_calc(listpath, basepath, outpath, title):
    """Execute everything."""
    print ('Calculating Similarity...')
    names = []
    scores = []
    authors = get_listlen(listpath)
    #Select one author
    for author in authors:
        row = []
        refauthor = get_authors(basepath,author)
        names.append(refauthor[0][0])
        #Compare all other authors to selected author
        for compare in authors:
            testauthor = get_authors(basepath, compare)
            vectorspace = construct_vectorspace(refauthor[1:], testauthor[1:])
            vector1 = populate_vectorspace(refauthor, vectorspace)
            vector2 = populate_vectorspace(testauthor, vectorspace)
            sim_score = calc_cos(vector1, vector2)
            row.append(sim_score)
        scores.append(row)
    print ('\nCalculation Completed For ' + str(len(names)) + ' PIs')
    plot_matrix(title, names, scores, outpath)
    return 0

#Execute
if __name__ == '__main__':
    cos_calc(LISTPATH, BASEPATH, OUTPATH, TITLE)

