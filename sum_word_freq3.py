#Sum Freq Appearing Words Across Multiple Papers
#Generate distribution of word frequency in pubmed abstracts
#Read in list of PIs as .txt
#Read in list of trivial words to exclude

#Mine pubmed abstracts using biopython package and Entrez

#Interpret and extract abstract text from xml using ElementTree

#by matt stetz 03/2018

###

#eliminate integer division
from __future__ import division
#use biopython (Bio) to access Entrez from NCBI
from Bio import Entrez
#use xml.etree / ElementTree to handle XML format
from xml.etree import cElementTree as ET

import pylab as plt
import numpy as np

###

#Provide email (required for Entrez)
Entrez.email="mstetz@pennmedicine.upenn.edu"

#Datapath to txt file with list of author names
INPATH='/Users/matthewstetz/Desktop/BMB/authors.txt'

#Datapath to txt file with list of trivial words to exclude
TRPATH='/Users/matthewstetz/Desktop/BMB/trivial.txt'

#Datapath to output directory
OUTPATH='/Users/matthewstetz/Desktop/BMB/ABS/'

#Maximum number of records to check
RECNUM=500

#Maximum number of keywords to compare
WORDNUM=10

#Minimum frequency of word occurence (2 is good, bigger numbers results in too much filtering)
FREQ=2

###

#query pubmed to get list of papers from a given author
#returns list of pubmed IDs
def query(name,number):
     handle = Entrez.esearch(db="pubmed",term=name,retmax=number)
     tree = ET.parse(handle)
     root = tree.getroot()
     #Get PubMed IDs for all papers
     papers=[]
     for Id in root.iter('Id'):
          papers.append(Id.text)
     return papers

#fetch abstract for a single paper
def get_abstract(paper):
     handle = Entrez.efetch(db="pubmed",id=paper,retmode='xml')
     tree = ET.parse(handle)
     root = tree.getroot()
     #account for papers that do not have abstracts
     abstract=[]
     for ab in root.iter('AbstractText'):
          try:
               abstract=(ab.text).split()
          except:
               print 'Abstract text not readable!'
               pass
     return abstract

#get rid of any capital letters 
def no_caps(text):
     new_text=[]
     for word in text:
          new_text.append(word.lower())
     return new_text

#get rid of any quotes
def no_quotes(text):
     new_text=[]
     for word in text:
          if word[0]=='"' and word[-1]!='"':
               new_text.append(word[1:])
          elif word[0]=='"' and word[-1]=='"':
               new_text.append(word[1:-1])
          elif word[0]!='"' and word[-1]=='"':
               new_text.append(word[0:-1])
          else:
               new_text.append(word)
     return new_text

#get rid of any other punctuation
#other punctionan only shows up at end of word
def no_punctuation(text):
     new_text=[]
     for word in text:
          if word[-1]=='.' or word[-1]==',':
               new_text.append(word[0:-1])
          else:
               new_text.append(word)
     return new_text

#get rid of any trivial words that
#are used frequently e.g. 'and', 'but', 'by', 'with', etc
def no_trivialwords(text,trivial):
     new_text=[]
     for word in text:
          if word not in trivial:
               new_text.append(word)
     return new_text

#count occurance of words in abstract
def count(text):
     count=[]
     final=[]
     c=1
     for x in range(len(text)):
          count.append(text[x])
          i=0
          while i < len(text):
               if text[x]==text[i]: 
                    c+=1
                    i+=1
               else:
                    i+=1
          count.append(c)
          final.append(count)
          count=[]
          c=0
          i=0
     #order list by frequency of occurcence
     final.sort(key=lambda x: int(x[1]))
     return final

def get_freqwords(text,freq):
     words=[]
     for x in range(len(text)):
          if int(text[x][1]) >= freq:
               words.append(text[x])
          else:
               pass
     #check for rare case where abstract does not repeat a single word
     if len(words)==0:
          words=text
     else:
          pass
     return words

def remove_duplicates(text):
     temp=[]
     final=[]
     for x in range(len(text)):
          if text[x][0] not in temp:
               temp.append(text[x][0])
               final.append(text[x])
          else:
               pass
     return final

#run all cleaning functions 
def cleanup(paper,trivial):
     abstract=no_trivialwords(no_punctuation(no_quotes(no_caps(paper))),trivial)
     return abstract

#this function is used to extract cleaned abstracts used for mining
def run(paper,trivial,freq):
     #get the abstract
     abstract=get_abstract(paper)
     #clean the abstract 
     abstract=cleanup(abstract,trivial)
     #count freq of occurrence for each word
     abstract=count(abstract)
     #filter out words with low freq of occurence
     abstract=get_freqwords(abstract,freq)
     #remove multiple spearate counts of same word
     abstract=remove_duplicates(abstract)
     return abstract    

#sum up all instances of a single word
def summ(master,abstracts):
     temp=[]
     for x in range(len(master)):
          temp.append(master[x][0])
     for y in range(len(abstracts)):
          for p in range(len(master)):
               if abstracts[y][0]==master[p][0]:
                    master[p][1]=master[p][1]+abstracts[y][1]
               elif abstracts[y][0] not in temp:
                    master.append(abstracts[y])
                    temp.append(abstracts[y][0])
               else:
                    pass
     return master

#prevent sporadic disconnects from interrupting script
def no_httperror(x,papers,trivial,freq,master):
     try:
          abstract=run(papers[x],trivial,freq)
          update=summ(master,abstract)
          master=update
     except:
          print 'bad connection: retrying'
          abstract=run(papers[x],trivial,freq)
          update=summ(master,abstract)
          master=update
     return master

#sort and take top number of key words
def cull(master,amount):
     culled=[]
     master.sort(key=lambda x: int(x[1]))
     cull=master[-int(amount):]
     return cull

#remove dummy name used for list operations for combining plural/singular versions of the same word
def remove_xxx(words):
     new=[]
     for x in range(len(words)):
          if words[x][0]!='xxx':
               new.append(words[x])
          else:
               pass
     return new
                   
#check for words that are counted separately for singular and plural versions
#combine singular and plural words
#i know this is not the most elegant code
def no_plural(words):
     new=[]
     temp=[]
     for x in range(len(words)):
          for y in range(len(words)):
               if (words[x][0]==words[x][0][0:-1] or words[x][0][0:-1]==words[y][0]) and (words[x][0]!='xxx'):
                    words[x][1]=words[x][1]+words[y][1]
                    new.append(words[x])
                    temp.append(words[x][0])
                    words[y][0]='xxx'
               else:
                    pass
          if (words[x][0] not in temp) and (words[x][0]!='xxx'):
               new.append(words[x])
          else:
               pass
     new=remove_xxx(new)
     temp=[]
     words=[]
     return new

#for the unusual situation where the first abstract is empty
def get_master(papers,trivial,freq):
     i=0
     while i < len(papers):
          master=run(papers[i],trivial,freq)
          if len(master)!=0:
               master=master
               break
          else:
               i+=1
     return master,i

#clean up and organize the master list of mined keywords 
def clean_master(master,number):
     master.sort(key=lambda x: int(x[1]))
     master=cull(master,number)
     master=no_plural(master)
     master.sort(key=lambda x: int(x[1]))
     master.reverse()
     return master

#mine abstracts from pubmed using entrez
def mine(papers,trivial,freq,number):
     master,index=get_master(papers,trivial,freq)
     for x in range(len(papers)):
          if x!=index:
               print x
               try:
                    abstract=run(papers[x],trivial,freq)
                    update=summ(master,abstract)
                    master=update
               except:
                    print 'bad connection: retrying'
                    master=no_httperror(x,papers,trivial,freq,master)
     master=clean_master(master,number)
     return master

#plot results as formatted bar graph (optional)
def plot(master,name,outpath):
     x_axis=np.arange(len(master))
     x_names=[master[x][0] for x in range(len(master))]
     y_vals=[master[x][1] for x in range(len(master))]
     plt.bar(x_axis,y_vals,align='center')
     plt.xticks(x_axis,x_names,rotation=45,fontsize=16)
     plt.yticks(fontsize=16)
     plt.tight_layout()
     plt.savefig(outpath+name+'.png')
     plt.clf()
     return 0

#save results to organized text file
def printout(master,author,outpath):
     with open(outpath,'w') as text:
          text.write(author+'\n')
          pad=max(len(row[0]) for row in master)+5
          for x in range(len(master)):
               text.write(master[x][0])
               text.write(str(master[x][1]).rjust(pad-len(master[x][0]))+'\n')
     return 0

#read input text file and extract relevant data
def read_names(inpath):
     names=[]
     with open(inpath,'r') as text:
          text=text.readlines()
          for name in text:
               names.append(name.strip('\n'))
     return names

#execute everything
def main(inpath,trpath,outpath,freq,recnum,wordnum):
     author=read_names(inpath)
     trivial=read_names(trpath)     
     for name in author:
          print '\n'+name
          outpathtxt=outpath+name+'.txt'
          papers=query(name,recnum)
          try:
               master=mine(papers,trivial,freq,wordnum)
               printout(master,name,outpathtxt)
               plot(master,name,outpath)
          except UnboundLocalError:
               print 'no abstracts found for '+name
               print 'check spelling of author name'
     return 0

if __name__=='__main__':         
     main(INPATH,TRPATH,OUTPATH,FREQ,RECNUM,WORDNUM)
