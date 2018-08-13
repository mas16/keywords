#Calculate cosine similarity for bag of words
#Analysis of PI abstracts mined from PubMed/MedLine
#By Matt Stetz 05/2018

###

###Import Libraries
import matplotlib
#Need to change matplotlib backend for plotting
matplotlib.use('Agg')
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

###User input below

#Working directory path
BASEPATH='/Users/matthewstetz/Desktop/BMB/Abs/'

#List of PIs path
LISTPATH='/Users/matthewstetz/Desktop/BMB/authors.txt'

#Path to where you want to save the output to
OUTPATH='/Users/matthewstetz/Desktop/BMB/bmb.png'

###No more user input needed

###Functions
#Get list of authors from text
def get_listlen(listpath):
     with open(listpath) as authors:
          data=authors.read()
          data=data.split('\n')
     return data

#Get author keywords from text
def get_authors(basepath,authors):
     authorpath=basepath+authors+'.txt'
     keywords=get_keywords(authorpath)
     return keywords

#Parse keywords
def get_keywords(datapath):
     with open(datapath) as keywords:
          wordlist=keywords.read()
          wordlist=wordlist.split()
          wordlist=organize_keywords(wordlist)
     return wordlist

#Clean keyword vectors     
def organize_keywords(keywords):
     temp=[]
     final=[]
     for x in range(len(keywords)):
          if x % 2 != 0:
               temp.append(keywords[x-1])
               if x > 1:
                    temp.append(float(keywords[x]))
               else:
                    temp.append(keywords[x])
               final.append(temp)
               temp=[]
          temp=[]
     return final

#Construct text vector space for two authors
def construct_vectorspace(v1,v2):
     t1=[]
     for x in range(len(v1)):
          if v1[x][0] not in t1:
               t1.append(v1[x][0])
     for x in range(len(v2)):
          if v2[x][0] not in t1:
               t1.append(v2[x][0])
     return t1

#Populate text vector space with keyword frequency for each author
def populate_vectorspace(v1,combined):
     t1=[]
     for x in range(len(combined)):
          i=0
          c=0
          while i < len(v1):
               if combined[x]==v1[i][0]:
                    t1.append(v1[i])
                    i+=1
                    c+=1
               else:
                    i+=1
          if c==0:
               t1.append([combined[x],0])
     return t1

#Calculate cosine similarity score for two authors
def calc_cos(v1,v2):
     temp_p,final=[],[]
     temp_a,temp_b=[],[]
     for x in range(len(v1)):
          temp_p.append(v1[x][1]*v2[x][1])
          temp_a.append(v1[x][1]**2)
          temp_b.append(v2[x][1]**2)
     sum_p=sum(temp_p)
     sum_a=sum(temp_a)
     sum_b=sum(temp_b)
     sum_a=np.sqrt(sum_a)
     sum_b=np.sqrt(sum_b)
     return sum_p/(sum_a*sum_b)

#Plot using seaborn heatmap
def plot_matrix(names,data,outpath):
     title='Similarity Scores for BMB Primary Faculty\n'
     ax=sns.heatmap(data,cmap='coolwarm')
     ax.set_xticklabels(names,rotation='vertical')
     ax.set_yticklabels(names,rotation='horizontal')
     plt.tick_params(direction='out')
     plt.title(title,fontweight='bold')
     plt.tight_layout()
     plt.savefig(outpath)
     print '\nCorrelation Matrix Saved'
     return 0

#Main function
def cos_calc(listpath,basepath,outpath):
     print 'Calculating Similarity...'
     temp,final,names=[],[],[]
     authors=get_listlen(listpath)
     for x in range(len(authors)):
          refauthor=get_authors(basepath,authors[x])
          names.append(refauthor[0][0])
          i=0
          while i < len(authors):
               testauthor=get_authors(basepath,authors[i])
               vectorspace=construct_vectorspace(refauthor[1:],testauthor[1:])
               vector1=populate_vectorspace(refauthor,vectorspace)
               vector2=populate_vectorspace(testauthor,vectorspace)
               sim_score=calc_cos(vector1,vector2)
               temp.append(sim_score)
               i+=1
          final.append(temp)
          temp=[]
     print '\nCalculation Completed For '+str(len(names))+' PIs'
     plot_matrix(names,final,outpath)
     return 0

#Execute
if __name__=='__main__':
     cos_calc(LISTPATH,BASEPATH,OUTPATH)
