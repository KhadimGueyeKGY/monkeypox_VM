# -*- coding: utf-8 -*-

"""
@author: kgy
"""


import os
import pandas as pd
from datetime import  datetime
import webbrowser 



os.system("cd packages/ena-mpx-prep/ ; python metadata_processing.py")
os.system("mkdir -p packages/ena-mpx-prep/data/seq")
os.system("mkdir -p packages/monkeypox/data")
matadata = pd.read_csv("packages/ena-mpx-prep/data/ENA_Search.tsv",sep="\t")
id = matadata.accession
contry = matadata.country
date = matadata.collection_date
date_public = matadata.first_public

# #**************************************** Downloading and heading Fasta ....
print('\n\n> Downloading Fastas ... [ '+str(datetime.now())+" ]")
#os.system("mkdir -p data")
for i in range(len(id)):
    os.system('curl "https://www.ebi.ac.uk/ena/browser/api/fasta/'+id[i]+'?download=true" --output packages/ena-mpx-prep/data/seq/'+id[i]+'.fasta')
os.system("cat packages/ena-mpx-prep/data/seq/* > packages/ena-mpx-prep/data/sequences.fasta")

print('\n\n> Heading Fastas ... [ '+str(datetime.now())+" ]")
Input = open("packages/ena-mpx-prep/data/sequences.fasta","r").read().split("\n")
Output = open("packages/monkeypox/data/sequences.fasta","w")
index = 0
for i in range(len(Input)):  
    if Input[i].find(">") != -1 :
        Output.write(">"+id[index]+"\n")
        index += 1
    else:
        Output.write(Input[i]+"\n")
Output.close()

# #**************************************** MetaData prep...

def contry_to_continent(ContryInput):
    ContinentContry = pd.read_csv("packages/continent_contry.txt", sep="\t")
    Continent = ContinentContry["Continent"]
    Contry = ContinentContry["Country"]
    for i in range(len(Contry)):
        if Contry[i].lower().find(ContryInput.lower()) != -1 :
            return Continent[i] , Contry[i]
    if i == len(Contry)-1 :
        return " " , " "
    if ContryInput.find(" ") != -1:
        CI = ContryInput.split(" ")
        for i in range(len(Contry)):
            for j in range(len(CI)):
                if Contry[i].lower().find(CI[j].lower()) == -1 :
                    break
            if j == len(CI)-1 :
                return Continent[i] , Contry[i]
        if i == len(Contry) -1 :
            return  " " , " "
        

print('\n\n> MetaData prep ... [ '+str(datetime.now())+" ]")
data_ref = pd.read_csv("packages/monkeypox/data_oj/metadata.tsv",sep="\t")
head_ref = [str(i) for i in data_ref]
first_line_ref = open("packages/monkeypox/data_oj/metadata.tsv","r").read().split("\n")[1].split("\t")
meatadata_output = open("packages/monkeypox/data/metadata.tsv","w")
line = "\t".join(head_ref)
meatadata_output.write(line+"\n")
for i in range(len(id)):
    first_line_ref[head_ref.index("accession")]= id[i]
    first_line_ref[head_ref.index("date")]= date[i]
    first_line_ref[head_ref.index("date_submitted")]= date_public[i]
    if contry[i].find(":") != -1 :
        ctt, ctry = contry_to_continent(contry[i].split(":")[0])
        if ctry != " ":
            ctry = contry[i].split(":")[0]
        first_line_ref[head_ref.index("country")]= ctry
        first_line_ref[head_ref.index("division")]= contry[i].split(":")[1]
        first_line_ref[head_ref.index("region")]= ctt
    else:
        ctt, ctry = contry_to_continent(contry[i])
        if ctry != " ":
            ctry = contry[i]
        first_line_ref[head_ref.index("country")]= ctry
        first_line_ref[head_ref.index("region")]= ctt
    
    line = "\t".join([str(e) for e in first_line_ref])
    meatadata_output.write(line+"\n")
meatadata_output.close()


#************************************** Run analysis pipeline ...

print('\n\n> Run analysis Nextstrain pipeline ... [ '+str(datetime.now())+" ]")

os.system("cd packages/monkeypox ; nextstrain build --docker --cpus 1 . --configfile config/config_mpxv.yaml")
os.system("cp packages/monkeypox/auspice/* packages/aupice_res/aupice/")
os.system("sh rename_mpn.sh")
# os.system("cd packages/monkeypox ; nextstrain build --docker --cpus 1 . --configfile config/config_hmpxv1.yaml") 


#************************************** Visualize results ....

#print('\n\n> Visualize results ... [ '+str(datetime.now())+" ]") 
#os.system("cd packages/monkeypox ; nextstrain view auspice/") 

#webbrowser.open_new_tab("http://127.0.0.1:4000/monkeypox/mpxv")








	
