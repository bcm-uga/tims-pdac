#cd Desktop/Post-doc/Luke/projects/thema_surv/real_data
#python3.7
#python3 (currently 3.13)

import os
os.getcwd()

import pandas as pd
from bs4 import BeautifulSoup as BS
from urllib.request import urlopen
import webbrowser

import gseapy as gp
from gseapy.plot import barplot, dotplot

# LOAD LIST OF GENES FROM R
listing0 = pd.read_csv('results/07_genes_hit.csv', keep_default_na=False).values.tolist()

listing = []
for i in range(len(listing0)):
    if ';' in listing0[i][0]:
        listing.extend(listing0[i][0].split(';'))
    elif listing0[i][0]!='':
        listing.append(listing0[i][0])
len(listing)

# SCRAPING
# LIST OF URLS
urls = []
for i in range(len(listing)):
 urls.append('https://www.ncbi.nlm.nih.gov/pubmed?term=%22Humans%22%5BMesh%5D+AND+(pancreatic+cancer%5BMeSH+Terms%5D)+AND+%22' + listing[i] + '%22')

# NB OF ITEMS FOR EACH ALIAS
items = []
for url in urls:
    text = urlopen(url)
    html = text.read()
    soup = BS(html, "html.parser")
    res = soup.find("meta", {"name": "log_resultcount"})
    if res is None:
        items.append(1)
    else:
        items.append(soup.find("meta", {"name": "log_resultcount"})["content"])

# Manually look at each result page
#useful_url = []
#useful_aliases = []

#for i in range(len(items)):
#    if items[i] != 0:
#        useful_url.append(urls[i])

#for i in range(len(items)):
#    if items[i] != 0:
#            useful_aliases.append(liste_listing[i])

#for url in useful_url:
#    webbrowser.open(url, new=0)

#c = open("20180830_total_Pubmed_scraping.csv", "w")
#for number in total:
#    c.write(str(number))
#    c.write('\n')
#for gene in useful_genes:
#    c.write(gene)
#    c.write('\n')
#c.close()

# Export results
pd.DataFrame(list(zip(listing,items))).to_csv('results/07_genes_hit_pubmed.csv')

# ENRICHMENT ANALYSIS (MAG : j'ai commenté ce paragraphe car j'avais des problèmes pour faire tourner gseapy)
#gene_listing = listing0
#for i in range(len(gene_listing)):
#    if ';' in gene_listing[i][0]:
#        gene_listing[i] = gene_listing[i][0].split(';')

#enr = [gp.enrichr(gene_list=genes,
           #      gene_sets='GO_Biological_Process_2021',
            #     organism='Human',
             #    outdir='') for genes in gene_listing if genes[0] != '']
#print(enr[0].results[['Term', 'Adjusted P-value', 'Overlap']].head(10))
#barplot(enr[0],title='G0_2021',)

#modules = pd.read_csv('results/07_genes_hit.csv', keep_default_na=False).iloc[:,0].tolist()
#modules = [i for (i, v) in zip(modules, [genes[0] != '' for genes in gene_listing]) if v]
#_ = [enrich.results.to_csv('results/07_genes_enrich_'+module+'.csv') for (enrich,module) in zip(enr,modules)]


# USE PROTEIN ATLAS TO DETERMINE THE CELLULAR LOCATION FOR EACH PROTEIN
# OPEN GENECARDS FOR INTERESTING DMR
#interesting_DMR = ['DMR'+str(i) for i in [2,22,34,68]]
#interesting_genes = pd.DataFrame(listing0).loc[pd.DataFrame(listing0)[0].isin(interesting_DMR)][1].tolist()
#urls_genecards = []
#urls_protein = []
#for i in range(len(interesting_genes)):
 #   if ';' in interesting_genes[i]:
  #      urls_genecards.extend(['https://www.genecards.org/cgi-bin/carddisp.pl?gene='+gene for gene in interesting_genes[i].split(';')])
   #     urls_protein.append(['https://www.proteinatlas.org/search/'+gene.lower() for gene in interesting_genes[i].split(';')])
    #else:
     #   urls_genecards.append('https://www.genecards.org/cgi-bin/carddisp.pl?gene='+interesting_genes[i])
      #  urls_protein.append('https://www.proteinatlas.org/search/'+interesting_genes[i].lower())

#for url in urls_genecards:
   # webbrowser.open(url, new=0)