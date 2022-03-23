import re, string, unicodedata
import inflect
from nltk import word_tokenize, sent_tokenize
from nltk.corpus import stopwords
from nltk.stem import LancasterStemmer, WordNetLemmatizer
import nltk.data
from geniatagger import GENIATagger

from Bio import Entrez
import xml.etree.ElementTree as ET
import os
import itertools
from lxml import etree
from bs4 import BeautifulSoup
import re

class TestPmcFilter:
    def cleanhtml(self, raw_html):
        cleanr = re.compile('<.*?>')
        cleantext = re.sub(cleanr, '', raw_html)
        return cleantext

    def load_forbidden_expressions(self):
        exps=['mutate gene', 'mutation analysis','tumor mutation','larrge-analysis','large-scale','exome','rna-seq','microarray','gene expression','gene expression pattern','mutational spectra', 'gene set enrichment analysis','next-generation sequencing','gene expression profile']

        return exps

    def load_interacting_verbs(self):
        verbs=['binding','subunit','unphosphorylate','phosphorylation','colocalization', 'phosphorylate', 'bind', 'coimmunopurify','co-purify','coimmunopurified', 'cofractionate', 'cofractionates','downregulate', 'upregulate','colocalize', 'connection','complex','intermediate','interaction','activate', 'associate','cleave','crosslink','immunoprecipitate','interact','recruit','impair','induce','mediate','regulate','receptor','ligand']
        # removed: decrease, bind, signaling
        return verbs

    def load_experimental_methods(self):
        methods=[]

        f=open("sparql_detection_methods.tsv", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            methods.append(l[1].split(";"))
        f.close()

        return methods

    def organize_sentences_with_bait(self, file_evaluation_pairs, folder, file_pairs): 
        #p=Nlp_pre_processing_steps()
        tokenizer = nltk.data.load('tokenizers/punkt/PY3/english.pickle')
        executable_path = os.path.join(".", "geniatagger-3.0.2", "geniatagger")
        tagger = GENIATagger(executable_path)

        #lt = Literature_hits(folder+file_pairs)
        lt="bait"
        
        methods=self.load_experimental_methods()
        verbs=self.load_interacting_verbs()
        stop_expressions=self.load_forbidden_expressions()

        parser = etree.XMLParser(recover=True)

        pairs=[]

        if(os.path.isdir(folder+"processed_sentences")):
            # Loading processed pairs
            for f in os.listdir(folder+"processed_sentences"):
                if(f.startswith("scs_")):
                    pr=f.replace("scs_","").replace(".tsv","").split("-")
                    pr=[ pr[0].lower(), pr[1].lower() ]
                    pairs.append(pr)
        else:
            os.system("mkdir "+folder+"processed_sentences")

        if(os.path.isdir(folder+"pmc_articles") or os.path.isdir(folder+"pubmed_articles")):
            count=0
            f=open(folder+file_evaluation_pairs,"r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(int(l[2])>0 or int(l[4])>0):
                    p1=l[0].split("|")[1].lower()
                    p2=l[1].split("|")[1].lower()
                    
                    pair=[p1, p2]
                    count+=1
                    if( not pair in pairs ):
                        p1=p1
                        p2=p2
                        print("#", count, " Evaluating ", p1, " and ", p2)

                        identifier_result=p1+"-"+p2

                        g=open(folder+"processed_sentences/scs_"+identifier_result+".tsv","w")
                        g.write("id_pmc\tsentence location\thighlighted phrase\tinteracting verbs\tnumber of interacting verbs\texperimental methods found\tnumber of experimental methods\tprotein entities found\tnumber of protein entities\n")
                        g.close()
                        
                        type_="pmc"
                        pmids=[]
                        if(l[5]!=""):
                            pmids=l[5].split(",")
                        for a in pmids:
                            self.process_paper(folder, a, type_, identifier_result, tokenizer, tagger, lt, methods, verbs, stop_expressions, parser, pair)

                        type_="pubmed"
                        pmids=[]
                        if(l[3]!=""):
                            pmids=l[3].split(",")
                        for a in pmids:
                            self.process_paper(folder, a, type_, identifier_result, tokenizer, tagger, lt, methods, verbs, stop_expressions, parser, pair)

                        if(os.path.getsize(folder+"processed_sentences/scs_"+identifier_result+".tsv")==198):
                            os.system("rm "+folder+"processed_sentences/scs_"+identifier_result+".tsv")
                            
            f.close()
        else:
            print("Error: You have to run step 2 first because the pmc/pubmed papers found in step 1 are not processed and stored")

    def process_paper(self, folder, a, type_, identifier_result, tokenizer, tagger, lt, methods, verbs, stop_expressions, parser, pair):
        if( os.path.isfile(folder+type_+'_articles/'+a+'.xml') and os.path.getsize(folder+type_+'_articles/'+a+'.xml')<=50000000 ):
            #print(p1, p2, 'pmc_articles/'+a+'.xml')
            try:
                # fixing xml file containing possible undesirable characters for xml
                root = etree.parse(folder+type_+'_articles/'+a+'.xml', parser=parser)
                fixed_xml=str(etree.tostring(root)).replace("b\'","").replace("'","")
                h=open(folder+type_+'_articles/'+a+'.xml',"w")
                h.write(fixed_xml)
                h.close()

                tree = ET.parse(folder+type_+'_articles/'+a+'.xml')
                root = tree.getroot()
                
                stemmer = LancasterStemmer()

                # Removing unnecessary tags and their contents
                txt=ET.tostring(root, encoding="unicode")
                txt=txt.replace("\\n","")
                soup = BeautifulSoup(txt, "xml")
                for tag in soup.find_all(['table','table-wrap','table-wrap-foot','supplementary-material','xref','sup','fig']):
                    tag.decompose()
                h=open(folder+type_+'_articles/'+a+'.xml',"w")
                h.write(str(soup.prettify()))
                h.close()
        
                sentence_groups={}

                sentences=[]
                for tg in root.findall('abstract'):
                    txt=ET.tostring(tg, encoding="unicode") 
                    txt=txt.replace("\\n","").replace("\\","")
                    temp=tokenizer.tokenize(self.cleanhtml(txt))
                    #print(temp)
                    for s in temp:
                        tokenized = nltk.word_tokenize(s.replace("\n",""))
                        #normalized = p.normalize(tokenized)
                        sentences.append(" ".join(tokenized).replace("/"," / ").replace("-"," - "))
                sentence_groups["abstract"]=sentences

                sentences=[]
                for scs in root.findall('sentences'):
                    for tg in scs.findall('textgroup'):
                        txt=ET.tostring(tg, encoding="unicode") 
                        txt=txt.replace("\\n","").replace("\\","")
                        temp=tokenizer.tokenize(self.cleanhtml(txt))
                        #print(temp)
                        for s in temp:
                            tokenized = nltk.word_tokenize(s.replace("\n",""))
                            #normalized = p.normalize(tokenized)
                            sentences.append(" ".join(tokenized).replace("/"," / ").replace("-"," - "))
                sentence_groups["body"]=sentences

                for idst in sentence_groups.keys():
                    #print("---Processing ", len(sentence_groups[idst]), "sentences in ", idst)
                    for s in sentence_groups[idst]:
                        #print(s)
                        count_forbidden_expressions=0
                        for exps in stop_expressions:
                            if(s.lower().find(exps)!=-1):
                                count_forbidden_expressions+=1
                        if(count_forbidden_expressions==0):
                            found_methods=[]
                            found_interact_verbs=[]
                            p1=pair[0]
                            p2=pair[1]
                            
                            c={}
                            c[pair[0]]=0
                            c[pair[1]]=0
                            
                            w=[]
                            pos=[]
                            entities=[]
                            
                            tokens=s.split(" ")
                            for word in tokens:
                                
                                word=word.lower()
                                condition=(word.find(p1)!=-1)
                                if(condition):
                                    c[p1]+=1
                                    if(not word in entities):
                                        entities.append(word)
                                    
                                condition=(word.find(p2)!=-1)
                                if(condition):
                                    c[p2]+=1
                                    if(not word in entities):
                                        entities.append(word)
                                
                                base_form=stemmer.stem(word)
                                if(base_form in verbs):
                                    found_interact_verbs.append(base_form)
                                
                                for met in methods:
                                    if(base_form in met):
                                        found_methods.append(base_form)
                                        break
                                    
                                w.append(word)
                                    
                            """for word, base_form, pos_tag, chunk, named_entity in tagger.tag(s):
                                if(named_entity.find("protein")!=-1):
                                    entities.append(word)

                                if(lt!=None):
                                    p1=pair[0]
                                    p2=pair[1]

                                    word=word.lower()
                                    condition=(word.find(p1)!=-1)
                                    #if(p1 in lt.symbols.keys()):
                                    #    condition=(word.find(p1)!=-1 or word in lt.symbols[p1])

                                    if(named_entity.find("protein")!=-1 and condition ):
                                        c1+=1
                                    
                                    condition=(word.find(p2)!=-1)
                                    #if(p2 in lt.symbols.keys()):
                                    #    condition=(word.find(p2)!=-1 or word in lt.symbols[p2])
                                        
                                    if(named_entity.find("protein")!=-1 and condition ):
                                        c2+=1
                                
                                if(base_form in verbs):
                                    found_interact_verbs.append(base_form)

                                if(named_entity.find("protein")==-1):
                                    for met in methods:
                                        if(base_form in met):
                                            found_methods.append(base_form)
                                            break
                                
                                w.append(base_form)
                                #pos.append(pos_tag)"""
                            
                            if(lt=="bait"):
                                main_condition = (c[p1]>0 and c[p2]>0)
                            else:
                                main_condition = (len(entities)>=2)
                                
                            if( main_condition  and len(found_interact_verbs)>0 ):
                                for word, base_form, pos_tag, chunk, named_entity in tagger.tag(s):
                                    if(named_entity.find("protein")!=-1 and (not word.lower() in entities)):
                                        entities.append(word.lower())
                                        
                                if(len(entities)<18 and (p1 in entities) and (p2 in entities)):
                                    with open(folder+"processed_sentences/scs_"+identifier_result+".tsv","a") as fg:
                                        #fg.write("%s\t%s\t%s\t%s\n" %(a, str(w), str(pos), str(entities) ) )
                                        fg.write("%s\t%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\n" %(a, idst, " ".join(w), "|".join(found_interact_verbs), len(found_interact_verbs), "|".join(found_methods), len(found_methods), "|".join(entities), len(entities) ) )
            except:
                print("not ok")
                
    def run(self, file_evaluation_pairs, folder, file_pairs):
        self.organize_sentences_with_bait(file_evaluation_pairs, folder, file_pairs)
        
import sys
folder="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/predppi_new_paper/ppi_pubminer_evaluation/"
fe="literature_evaluation_pairs.tsv"
fp="filtered_pairs_by_assocrules.tsv"
a=TestPmcFilter()
a.run(fe, folder, fp)