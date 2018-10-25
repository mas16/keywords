"""
Scrape NCBI database using Entrez
and generate bag of words model from xml.

Calculate TF and determine cosine similarity
to analyze research content networks.

By MAS 03/2018.

Updated by MAS 10/2018 - Style clean up.
"""

###
#Import libraries

#Eliminate integer division
from __future__ import division
#Use biopython (Bio) to use Entrez for scraping NCBI
from Bio import Entrez
#Use xml.etree / ElementTree to handle XML format
from xml.etree import cElementTree as ET
#For error handling during scraping
from urllib2 import HTTPError, URLError

import matplotlib
#Change matplotlib backend for plotting to avoid warning
matplotlib.use('Agg')

import pylab as plt
import numpy as np

###
#Define constants

#Provide email (required for Entrez)
Entrez.email = 'mstetz@pennmedicine.upenn.edu'

#Datapath to .txt file with list of author names
INPATH = '/Users/matthewstetz/Desktop/datascience/JHU_TEST/authorsbarrick.txt'

#Datapath to .txt file with list of trivial words to exclude
TRPATH = '/Users/matthewstetz/Desktop/datascience/JHU_TEST/trivial.txt'

#Datapath to output directory
OUTPATH = '/Users/matthewstetz/Desktop/datascience/JHU_TEST/'

#Optional: University Affiliation (may be necessary for common names)
AFFILIATION = ''

#Maximum number of records to check
RECNUM = 500

#Maximum number of keywords to compare
WORDNUM = 10

#Minimum frequency of word occurence
#(2 is good, bigger numbers results in too much filtering)
FREQ = 2

###
#Define functions

def query(name, number):
    """Query pubmed to get all paper IDs for a given author."""
    handle = Entrez.esearch(db='pubmed', term=name, retmax=number)
    tree = ET.parse(handle)
    root = tree.getroot()
    papers = [Id.text for Id in root.iter('Id')]
    return papers

def get_abstract(paper):
    """Fetch abstract for each paper ID."""
    handle = Entrez.efetch(db='pubmed', id=paper, retmode='xml')
    tree = ET.parse(handle)
    root = tree.getroot()
    """Code below accounts for papers that do not have an abstract."""
    abstract = []
    for ab in root.iter('AbstractText'):
        try:
            abstract = (ab.text).split()[:]
        except:
            print ('No abstract found for record!')
            raise
    return abstract

def no_caps(text):
    """Remove capital letters."""
    new_text = [word.lower() for word in text]
    return new_text

def no_quotes(text):
    """Remove quotation marks."""
    new_text = []
    for word in text:
        if word[0] == '"' and word[-1] != '"':
            new_text.append(word[1:])
        elif word[0] == '"' and word[-1] == '"':
            new_text.append(word[1:-1])
        elif word[0] != '"' and word[-1] == '"':
            new_text.append(word[0:-1])
        else:
            new_text.append(word)
    return new_text

def no_punctuation(text):
    """Remove punctuation at end of word."""
    new_text = []
    for word in text:
        if word[-1] == '.' or word[-1] == ',':
            new_text.append(word[0:-1])
        else:
            new_text.append(word)
    return new_text

def no_trivialwords(text, trivial):
    """Remove common trivial words (stopwords in sci-kit learn)."""
    new_text = [word for word in text if word not in trivial]
    return new_text

def count(text):
    """Determine TFs."""
    counts = []
    for entry in text:
        c = 0
        for subentry in text:
            if subentry == entry:
                c += 1
        counts.append([entry,c])
    counts.sort(key=lambda x: int(x[1]))
    return counts

def get_freqwords(text, freq):
    """Filter words by threshold TF."""
    words = [entry for entry in text if entry[1] >= freq]
    """Check for rare case where all TFs = 1."""
    if len(words) == 0:
        words = text
    return words
 
def remove_duplicates(text):
    """Remove duplicate occurences"""
    repo = []
    final = []
    for entry in text:
        if entry[0] not in repo:
            repo.append(entry[0])
            final.append(entry)
    return final

#Run all cleaning functions 
def cleanup(paper, trivial):
    abstract = no_caps(paper)
    abstract = no_quotes(abstract)[:]
    abstract = no_punctuation(abstract)[:]
    abstract = no_trivialwords(abstract, trivial)[:]
    return abstract

#This function is used to extract cleaned abstracts used for mining
def run(paper, trivial, freq):
    abstract = get_abstract(paper)
    abstract = cleanup(abstract, trivial)[:]
    abstract = count(abstract)[:]
    abstract = get_freqwords(abstract, freq)[:]
    abstract = remove_duplicates(abstract)[:]
    return abstract    

def summ(master, abstracts):
    """Get total TF across all abstracts for a given word."""
    ref_list = [entry[0] for entry in master]
    """Sum all instances of the word."""
    for entry in abstracts:
        for p in range(len(master)):
            if entry[0] == master[p][0]:
                master[p][1] = master[p][1] + entry[1]
            elif entry[0] not in ref_list:
                master.append(entry)
                ref_list.append(entry[0])
    return master

def no_httperror(x, papers, trivial, freq, master):
    """Prevent sporadic internet disconnects from interrupting script"""
    try:
        abstract = run(papers[x], trivial, freq)
        update = summ(master, abstract)
        master = update[:]
    except URLError:
        print ('Bad Internet Connection: Retrying')
        master = no_httperror(x, papers, trivial, freq, master)[:]
    except HTTPError:
        print ('HTTP Error: Retrying') 
        master = no_httperror(x, papers, trivial, freq, master)[:]
    return master

def cull(master, amount):
    """Sort and take top number of key words."""
    master.sort(key=lambda x: int(x[1]))
    cull = master[-int(amount):]
    return cull

def remove_xxx(words):
    """Remove dummy name used for list operations -
    for combining plural/singular versions of the same word.
    """
    new = [word for word in words if word[0] != 'xxx']
    return new

def no_plural(words):
    """Check for singular and plural versions of same word -
    If found, combine TFs.
    """
    new = []
    temp = []
    for word in words:
        for subword in words:
            if (word[0] == word[0][0:-1] or word[0][0:-1] == subword[0]) and (word[0] != 'xxx'):
                word[1] = word[1]+subword[1]
                new.append(word)
                temp.append(word[0])
                subword[0] = 'xxx'
        if (word[0] not in temp) and (word[0] != 'xxx'):
            new.append(word)
    new = remove_xxx(new)
    return new

def get_master(papers, trivial, freq):
    """Get a paper to start the TF count."""
    for x, paper in enumerate(papers):
        master = run(paper, trivial, freq)
        if len(master) != 0:
            master = master[:]
            break
    return master, x

#Clean up and organize the master list of mined keywords 
def clean_master(master, number):
    master.sort(key=lambda x: int(x[1]))
    master = cull(master,number)[:]
    master = no_plural(master)[:]
    master.sort(key=lambda x: int(x[1]))
    master.reverse()
    return master

def mine(papers, trivial, freq, number):
    """Scrape NCBI database."""
    master, index = get_master(papers, trivial, freq)
    print ('Scraping ...')
    for x, paper in enumerate(papers):
        if x != index:
            """Make sure the script can continue to run if
            connection gets interupted.
            """
            master = no_httperror(x, papers, trivial, freq, master)[:]
            if x%10 == 0:
                print ('...')
    master = clean_master(master, number)[:]
    return master

def plot(master, name, outpath):
    """Plot TFs as formatted bar plot."""
    x_axis = np.arange(len(master))
    x_names = [entry[0] for entry in master]
    y_vals = [entry[1] for entry in master]
    plt.bar(x_axis, y_vals, align='center')
    plt.xticks(x_axis, x_names, rotation=45, fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig(outpath+name+'.png')
    plt.clf()
    return 0

def printout(master,author,outpath):
    """Save results as .txt file."""
    with open(outpath,'w') as text:
        text.write(author+'\n')
        pad = max(len(row[0]) for row in master)+5
        for entry in master:
            try:
                text.write(entry[0])
                text.write(str(entry[1]).rjust(pad-len(entry[0]))+'\n')
            except UnicodeEncodeError:
                print ('Unicode Encode Error: Skipping')
                pass
    return 0

def read_names(inpath):
    """Read input data and extract names."""
    with open(inpath,'r') as text:
        text = text.readlines()
        names = [name.strip('\n') for name in text]
    return names

def main(inpath, trpath, outpath, freq, recnum, wordnum, affiliation):
    """Execute everything."""
    author = read_names(inpath)
    trivial = read_names(trpath)     

    for name in author:
        print ('\n'+name)
        
        outpathtxt = outpath+name+'.txt'

        if len(affiliation) > 0:
            search_name = name+' AND '+affiliation
            print (search_name)
        else:
            search_name = name
            
        try:
            papers = query(search_name, recnum)
        except URLError:
            print ('\n**ERROR: Bad Internet Connection**\n')
            
        try:
            master = mine(papers, trivial, freq, wordnum)
            printout(master, name, outpathtxt)
            plot(master, name, outpath)
            print ('Plotting TFs...')
        except UnboundLocalError:
            print ('Warning: No abstracts found for '+name)
            
    print ('\nScript Finished')
    return 0

if __name__ == '__main__':         
    main(INPATH, TRPATH, OUTPATH, FREQ, RECNUM, WORDNUM, AFFILIATION)
