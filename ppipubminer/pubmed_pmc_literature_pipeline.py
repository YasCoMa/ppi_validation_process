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
class Literature_hits:

    def __init__(self, pairsfile=None):
        if(pairsfile!=None):
            self.symbols=self.build_dict_other_symbols(pairsfile)

    def build_dictionaries_mapping(self):
        hgnc={}
        f=open("mapping_hgnc_uniprot.txt","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            hgnc[l[1]] = l[0]
        f.close()

        return hgnc

    def build_dict_other_symbols(self, pairsfile):
        mapp={}

        hgnc=self.build_dictionaries_mapping()

        proteins=[]
        f=open(pairsfile,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[0] in hgnc.keys() and l[1] in hgnc.keys()):
                if(not hgnc[l[0]] in proteins):
                    proteins.append(hgnc[l[0]].lower())

                if(not hgnc[l[1]] in proteins):
                    proteins.append(hgnc[l[1]].lower())
                
        f.close()

        print("Building symbols...")
        f=open("genenames_data.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[1] in proteins):
                prots=[l[1]]
                if(l[2]!=""):
                    prots += l[2].lower().split(", ")
                if(l[3]!=""):
                    prots += l[3].lower().split(", ")

                variations=[]
                for p in prots:
                    variations.append(''.join(c[0] for c in itertools.groupby(p)))
                    variations.append("p"+p)
                prots+=variations
                #print(prots)
                mapp[l[1]]=prots
        f.close()

        return mapp

    def get_pmcHits_for_pair(self, pair):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        
        if(pair[0] in self.symbols.keys() and pair[1] in self.symbols.keys()):
            names=[ "("+(" or ".join(self.symbols[pair[0]]))+")", "("+(" or ".join(self.symbols[pair[1]]))+")" ]
        else:
            names=pair
         
        handle = Entrez.esearch(db='pmc', sort='relevance', term="interaction ("+(" and ".join(names))+")" )
        results = Entrez.read(handle)
        hits = results['IdList']
        return hits

    # Store the title, abstract and sentences for each article (pmcid1.xml, ...) 
    def get_hits_on_pmc_sentences(self, hits, folder):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        new_hits=0
        
        id=",".join(hits)
        fetch = Entrez.efetch(db='pmc',resetmode='xml',id=id,rettype='full')
        with open(folder+'tempFile.xml', 'w') as f:
            f.write(fetch.read())
        tree = ET.parse(folder+'tempFile.xml')
        root = tree.getroot()
        
        if( not os.path.isdir(folder+"pmc_articles") ):
            os.system("mkdir "+folder+"pmc_articles")
        
        # article file creation
        newfile=[]
        pmids_order=[]
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for pmids in meta.findall('article-id'):
                        if(pmids.attrib['pub-id-type']=='pmc'):
                            if(not os.path.isfile(folder+"pmc_articles/"+pmids.text+".xml")):
                                newfile.append(True)
                                pmids_order.append(folder+"pmc_articles/"+pmids.text+".xml")
                                with open(folder+"pmc_articles/"+pmids.text+".xml", "a") as g:
                                    g.write("<article>")
                            else:
                                newfile.append(False)
                                pmids_order.append(folder+"pmc_articles/"+pmids.text+".xml")
                        
        #print(pmids_order)
        #print(len(pmids_order), len(newfile))
        # getting title
        c=0
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for group in meta.findall('title-group'):
                        for title in group.findall('article-title'):
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("<title>%s</title>" %(title.text) )
            c+=1

        # getting abstract
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<abstract>")
            
            c+=1
        c=0
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for abs in meta.findall('abstract'):
                        for sec in abs.findall('sec'):
                            for blocks2 in sec.findall('p'):
                                complete_text=ET.tostring(blocks2, encoding="unicode") 
                                if(newfile[c]):
                                    with open(pmids_order[c], "a") as g:
                                        g.write("%s" %(complete_text) )

                        for blocks in abs.findall('p'):
                            complete_text=ET.tostring(blocks, encoding="unicode") 
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("%s" %(complete_text) )
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</abstract>")
            c+=1

        # getting sentences
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<sentences>")
            
            c+=1
        c=0
        for article in root.findall('article'):
            for body in article.findall('body'):
                for section in body.findall('sec'):
                    complete_text=ET.tostring(section, encoding="unicode") 
                    if(newfile[c]):
                        with open(pmids_order[c], "a") as g:
                            g.write("<textgroup >%s</textgroup>" %( complete_text) )

                for blocks in body.findall('p'):
                    complete_text=ET.tostring(blocks, encoding="unicode") 
                    if(newfile[c]):
                        with open(pmids_order[c], "a") as g:
                            g.write("<textgroup >%s</textgroup>" %(complete_text) )

                    """for section2 in section.findall('sec'):
                        section_type2 = ""
                        if('sec-type' in section2.attrib):
                            section_type2 = "gtype='"+section2.attrib['sec-type']+"'"
                        for blocks2 in section2.findall('p'):    
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("<textgroup %s >%s</textgroup>" %(section_type2, blocks2.text) )

                    section_type = ""
                    if('sec-type' in section.attrib):
                        section_type = "gtype='"+section.attrib['sec-type']+"'"
                    for blocks in section.findall('p'):
                        if(newfile[c]):
                            with open(pmids_order[c], "a") as g:
                                g.write("<textgroup %s >%s</textgroup>" %(section_type, blocks.text) )"""
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</sentences>")
            c+=1

        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</article>")
            c+=1

        return new_hits

    def get_pubmedHits_for_pair(self, pair):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        
        if(pair[0] in self.symbols.keys() and pair[1] in self.symbols.keys()):
            names=[ "("+(" or ".join(self.symbols[pair[0]]))+")", "("+(" or ".join(self.symbols[pair[1]]))+")" ]
        else:
            names=pair
         
        handle = Entrez.esearch(db='pubmed', sort='relevance', term="interaction ("+(" and ".join(names))+")" )
        results = Entrez.read(handle)
        hits = results['IdList']
        return hits

    # Store the title, abstract and sentences for each article (pmcid1.xml, ...) 
    def get_hits_on_pubmed_sentences(self, hits, folder):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        new_hits=0
        
        id=",".join(hits)
        fetch = Entrez.efetch(db='pubmed',resetmode='xml',id=id,rettype='full')
        with open(folder+'tempFile.xml', 'w') as f:
            f.write(fetch.read())
        tree = ET.parse(folder+'tempFile.xml')
        root = tree.getroot()
        
        if( not os.path.isdir(folder+"pubmed_articles") ):
            os.system("mkdir "+folder+"pubmed_articles")
        
        # article file creation
        newfile=[]
        pmids_order=[]
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for pmids in front.findall('PMID'):
                    if(not os.path.isfile(folder+"pubmed_articles/"+pmids.text+".xml")):
                        newfile.append(True)
                        pmids_order.append(folder+"pubmed_articles/"+pmids.text+".xml")
                        with open(folder+"pubmed_articles/"+pmids.text+".xml", "a") as g:
                            g.write("<article>")
                    else:
                        newfile.append(False)
                        pmids_order.append(folder+"pubmed_articles/"+pmids.text+".xml")
                        
        #print(pmids_order)
        #print(len(pmids_order), len(newfile))
        # getting title
        c=0
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for meta in front.findall('Article'):
                    for title in meta.findall('ArticleTitle'):
                        if(newfile[c]):
                            with open(pmids_order[c], "a") as g:
                                g.write("<title>%s</title>" %(title.text) )
            c+=1

        # getting abstract
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<abstract>")
            
            c+=1
        c=0
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for meta in front.findall('Article'):
                    for abs_ in meta.findall('Abstract'):
                        for sec in abs_.findall('AbstractText'):
                            complete_text=ET.tostring(sec, encoding="unicode") 
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("%s" %(complete_text) )
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</abstract>")
            c+=1

        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</article>")
            c+=1

        return new_hits

    # evaluate sentences in the text-groups and save those that have the couple of proteins under evaluation

class Nlp_pre_processing_steps:
	## BEGIN NLP pre processing steps
	def remove_non_ascii(self, words):
		"""Remove non-ASCII characters from list of tokenized words"""
		new_words = []
		for word in words:
			new_word = unicodedata.normalize('NFKD', word).encode('ascii', 'ignore').decode('utf-8', 'ignore')
			new_words.append(new_word)
		return new_words

	def to_lowercase(self, words):
		"""Convert all characters to lowercase from list of tokenized words"""
		new_words = []
		for word in words:
			new_word = word.lower()
			new_words.append(new_word)
		return new_words

	def remove_punctuation(self, words):
		"""Remove punctuation from list of tokenized words"""
		new_words = []
		for word in words:
			new_word = re.sub(r'[^\w\s]', '', word)
			if new_word != '':
				new_words.append(new_word)
		return new_words

	def replace_numbers(self, words):
		"""Replace all interger occurrences in list of tokenized words with textual representation"""
		p = inflect.engine()
		new_words = []
		for word in words:
			if word.isdigit():
				new_word = p.number_to_words(word)
				new_words.append(new_word)
			else:
				new_words.append(word)
		return new_words

	def remove_stopwords(self, words):
		"""Remove stop words from list of tokenized words"""
		new_words = []
		for word in words:
			if word not in stopwords.words('english'):
				new_words.append(word)
		return new_words

	def stem_words(self, words):
		"""Stem words in list of tokenized words"""
		stemmer = LancasterStemmer()
		stems = []
		for word in words:
			stem = stemmer.stem(word)
			stems.append(stem)
		return stems

	def lemmatize_verbs(self, words):
		"""Lemmatize verbs in list of tokenized words"""
		lemmatizer = WordNetLemmatizer()
		lemmas = []
		for word in words:
			lemma = lemmatizer.lemmatize(word, pos='v')
			lemmas.append(lemma)
		return lemmas

	def normalize(self, words):
		words = self.remove_non_ascii(words)
		words = self.to_lowercase(words)
		words = self.remove_punctuation(words)
		words = self.replace_numbers(words)
		words = self.remove_stopwords(words)
		return words

	## END NLP pre processing steps

from lxml import etree
from bs4 import BeautifulSoup
import re
class Text_processing:
    def cleanhtml(self, raw_html):
        cleanr = re.compile('<.*?>')
        cleantext = re.sub(cleanr, '', raw_html)
        return cleantext

    def load_forbidden_expressions(self):
        exps=['mutate gene', 'mutation analysis','tumor mutation','larrge-analysis','larrge-scale','exome','rna-seq','microarray','gene expression','gene expression pattern','mutational spectra', 'gene set enrichment analysis','next-generation sequencing','gene expression profile']

        return exps

    def load_interacting_verbs(self):
        verbs=['binding','subunit','unphosphorylate','phosphorylation','colocalization', 'phosphorylate', 'bind', 'coimmunopurify','co-purify','coimmunopurified', 'cofractionate', 'cofractionates','downregulate', 'upregulate','colocalize', 'connection','complex','intermediate','interaction','activate', 'associate','cleave','crosslink','immunoprecipitate','interact','recruit','block','impair','induce','mediate','regulate']
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

    def organize_sentences_without_bait(self, folder, file_articles): # articles list has two columns, one with type (pmc or pubmed) and the other with article id
        tokenizer = nltk.data.load('tokenizers/punkt/PY3/english.pickle')
        executable_path = os.path.join(".", "geniatagger-3.0.2", "geniatagger")
        tagger = GENIATagger(executable_path)

        methods=self.load_experimental_methods()
        verbs=self.load_interacting_verbs()
        stop_expressions=self.load_forbidden_expressions()

        parser = etree.XMLParser(recover=True)

        if(not os.path.isdir(folder+"processed_sentences")):
            os.system("mkdir "+folder+"processed_sentences")
        
        identifier_result="general_mode"

        g=open(folder+"processed_sentences/scs_"+identifier_result+".tsv","w")
        g.write("id_pmc\tsentence location\thighlighted phrase\tinteracting verbs\tnumber of interacting verbs\texperimental methods found\tnumber of experimental methods\tprotein entities found\tnumber of protein entities\n")
        g.close()

        if(os.path.isdir(folder+"pmc_articles") or os.path.isdir(folder+"pubmed_articles")): 
            f=open(folder+file_articles,"r")
            for line in f:
                l=line.replace("\n","").split("\t")
                type_=l[0]
                a=l[1]
                lt=None
                pair=None

                self.process_paper(folder, a, type_, identifier_result, tokenizer, tagger, lt, methods, verbs, stop_expressions, parser, pair)
                            
            f.close()
        else:
            print("Error: You have to run step 1 first because the pmc/pubmed papers are not processed and stored")

    def organize_sentences_with_bait(self, file_evaluation_pairs, folder, file_pairs): 
        p=Nlp_pre_processing_steps()
        tokenizer = nltk.data.load('tokenizers/punkt/PY3/english.pickle')
        executable_path = os.path.join(".", "geniatagger-3.0.2", "geniatagger")
        tagger = GENIATagger(executable_path)

        lt = Literature_hits(folder+file_pairs)

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
                    pr=[ pr[0], pr[1] ]
                    pairs.append(pr)
        else:
            os.system("mkdir "+folder+"processed_sentences")

        if(os.path.isdir(folder+"pmc_articles") or os.path.isdir(folder+"pubmed_articles")):
            f=open(folder+file_evaluation_pairs,"r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(int(l[2])>0 or int(l[4])>0):
                    p1=l[0].split("|")[1]
                    p2=l[1].split("|")[1]
                    
                    pair=[p1, p2]
                    if( not pair in pairs ):
                        p1=p1.lower()
                        p2=p2.lower()

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

                        if(os.path.getsize(folder+"processed_sentences/scs_"+identifier_result+".tsv")==130):
                            os.system("rm "+folder+"processed_sentences/scs_"+identifier_result+".tsv")
                            
            f.close()
        else:
            print("Error: You have to run step 2 first because the pmc/pubmed papers found in step 1 are not processed and stored")

    def process_paper(self, folder, a, type_, identifier_result, tokenizer, tagger, lt, methods, verbs, stop_expressions, parser, pair):
        if( os.path.isfile(folder+type_+'_articles/'+a+'.xml') ):
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
                    for s in sentence_groups[idst]:
                        #print(s)
                        count_forbidden_expressions=0
                        for exps in stop_expressions:
                            if(s.lower().find(exps)!=-1):
                                count_forbidden_expressions+=1
                        if(count_forbidden_expressions==0):
                            found_methods=[]
                            found_interact_verbs=[]
                            c1=0
                            c2=0
                            
                            w=[]
                            pos=[]
                            entities=[]
                            for word, base_form, pos_tag, chunk, named_entity in tagger.tag(s):
                                if(named_entity.find("protein")!=-1):
                                    entities.append(word)

                                if(lt!=None):
                                    p1=pair[0]
                                    p2=pair[1]

                                    word=word.lower()
                                    condition=(word.find(p1)!=-1)
                                    if(p1 in lt.symbols.keys()):
                                        condition=(word.find(p1)!=-1 or word in lt.symbols[p1])

                                    if(named_entity.find("protein")!=-1 and condition ):
                                        c1+=1
                                    
                                    condition=(word.find(p2)!=-1)
                                    if(p2 in lt.symbols.keys()):
                                        condition=(word.find(p2)!=-1 or word in lt.symbols[p2])
                                        
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
                                #pos.append(pos_tag)
                            
                            if(lt!=None):
                                main_condition = (c1>0 and c2>0)
                            else:
                                main_condition = (len(entities)>=2)

                            if( main_condition and len(entities)<15 and len(found_interact_verbs)>0):
                                with open(folder+"processed_sentences/scs_"+identifier_result+".tsv","a") as fg:
                                    #fg.write("%s\t%s\t%s\t%s\n" %(a, str(w), str(pos), str(entities) ) )
                                    fg.write("%s\t%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\n" %(a, idst, "|".join(w), "|".join(found_interact_verbs), len(found_interact_verbs), "|".join(found_methods), len(found_methods), "|".join(entities), len(entities) ) )
            except:
                print("not ok")

class Running_config:

    def run_step1_mode1(self, file_pairs, folder):
        print("Step 1 - Initializing aliases of protein symbols...")
        lt=Literature_hits(folder+file_pairs)

        hgnc=lt.build_dictionaries_mapping()

        if( os.path.isdir(folder+"pmc_articles") ):
            os.system("rm "+folder+"pmc_articles/*")

        f=open(folder+"literature_evaluation_pairs.tsv","w")
        f.close()
        c=0
        f=open(folder+file_pairs,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[0] in hgnc.keys() and l[1] in hgnc.keys()):
                c+=1
                hits_pmc=lt.get_pmcHits_for_pair([hgnc[l[0]], hgnc[l[1]]])
                hits_pubmed=lt.get_pubmedHits_for_pair([hgnc[l[0]], hgnc[l[1]]])
                with open(folder+"literature_evaluation_pairs.tsv","a") as g:
                    g.write("%s\t%s\t%i\t%s\t%i\t%s\n" %(l[0]+"|"+hgnc[l[0]], l[1]+"|"+hgnc[l[1]], len(hits_pubmed), (",".join(hits_pubmed)), len(hits_pmc), (",".join(hits_pmc)) ) )

        f.close()

    def run_step2_mode1(self, file_evaluation_pairs, folder, file_pairs):
        print("Step 2 - Initializing aliases of protein symbols...")
        lit=Literature_hits(folder+file_pairs)
        
        print("Step 2 - Getting articles content...")
        f=open(folder+file_evaluation_pairs,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if( int(l[4])>0):
                hits=l[5].split(",")
                hits_sentences=lit.get_hits_on_pmc_sentences(hits, folder)

            if( int(l[2])>0):
                hits=l[3].split(",")
                hits_sentences=lit.get_hits_on_pubmed_sentences(hits, folder)

        f.close()

    def run_step3_mode1(self, file_evaluation_pairs, folder, file_pairs):
        t=Text_processing()
        t.organize_sentences_with_bait(file_evaluation_pairs, folder, file_pairs)

    def run_step1_mode2(self, folder, file_articles):
        lit=Literature_hits(None)
        
        f=open(folder+file_articles,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            hits=[l[1]]
            if(l[0]=="pmc"):
                hits_sentences=lit.get_hits_on_pmc_sentences(hits, folder)

            if(l[0]=="pubmed"):
                hits_sentences=lit.get_hits_on_pubmed_sentences(hits, folder)

        f.close()

    def run_step2_mode2(self, folder, file_articles):
        t=Text_processing()
        t.organize_sentences_without_bait(folder, file_articles)

    def run(self, args):
        if(args.folder!="" and os.path.isdir(args.folder)):
            if(args.execution_mode in [1,2]):
                if(args.execution_mode==1):
                    run=0
                    if(args.running_type_mode_1=="" ):
                        run=0
                    else:
                        if(args.running_type_mode_1 in [0,1,2,3]):
                            run=args.running_type_mode_1
                        else:
                            print("Error: invalid choice")

                    if(run==0 or run==1):
                        if(args.file_pairs==""):
                            print("Error: you have to give the file with pairs")
                        else:
                            print("Running step 1")
                            self.run_step1_mode1(args.file_pairs, args.folder)
                            if(run==0):
                                print("Running step 2")
                                self.run_step2_mode1("literature_evaluation_pairs.tsv", args.folder, args.file_pairs)
                                print("Running step 3")
                                self.run_step3_mode1("literature_evaluation_pairs.tsv", args.folder, args.file_pairs)
                    else:
                        if(args.file_evaluation==""):
                            print("Error: you have to give the evaluation file exported in step 1")
                        else:
                            if(args.file_evaluation==""):
                                print("Error: you have to give the file with pairs")
                            else:
                                if(run==2):
                                    self.run_step2_mode1(args.file_evaluation, args.folder, args.file_pairs)
                                else:
                                    self.run_step3_mode1(args.file_evaluation, args.folder, args.file_pairs)

                if(args.execution_mode==2):
                    run=0
                    if(args.running_type_mode_2=="" ):
                        run=0
                    else:
                        if(args.running_type_mode_2 in [0,1,2]):
                            run=args.running_type_mode_2
                        else:
                            print("Error: invalid choice")

                    if(run==0 or run==1):
                        if(args.file_articles==""):
                            print("Error: you have to give the file with articles")
                        else:
                            print("Running step 1")
                            self.run_step1_mode2(args.folder, args.file_articles)
                            if(run==0):
                                print("Running step 2")
                                self.run_step2_mode2(args.folder, args.file_articles)
                    else:
                        self.run_step2_mode2(args.folder, args.file_articles)
                
            else:
                print("Error: invalid choice of execution mode")
        else:
            print("Error: You have to specify a valid folder to store files")


import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' PPIPubMiner - Literature pipeline to find protein interactions on pubmed articles', formatter_class=RawTextHelpFormatter)
parser.add_argument("-em", "--execution_mode", action="store", help="1 - Mode using a list of protein pairs as bait\n\
2 - Mode that tries to find sentences of PPI context for any protein pairs given a list of articles", type=int)
parser.add_argument("-fo", "--folder", action="store", help="(For both modes) Folder to store the files (use the folder where the required files can be found, ex.: /home/user/experiment/ )\n")

parser.add_argument("-rtm1", "--running_type_mode_1", action="store", help="(For mode 1) 0 (default) - Run all steps\n\
1 - Run step 1 (Get mentions of both proteins in PMC articles)\n\
2 - Run step 2 (Get the PMC or Pubmed files, clean and store them)\n\
3 - Run step 3 (Get the exact sentences where the proteins were found on interacting context)", type=int)
parser.add_argument("-fp", "--file_pairs", action="store", help="(For mode 1) File with the pairs (two columns with uniprot identifiers in tsv format)")
parser.add_argument("-fe", "--file_evaluation", action="store", help="(For mode 1) File exported after step 1 execution in tsv format\n")

parser.add_argument("-rtm2", "--running_type_mode_2", action="store", help="(For mode 2) 0 (default) - Run all steps\n\
1 - Run step 1 (Get the PMC or Pubmed files from the given list, clean and store them)\n\
2 - Run step 2 (Get the exact sentences where the proteins were found on an interacting context)", type=int)
parser.add_argument("-fa", "--file_articles", action="store", help="(For mode 2) File with the articles (First column indicating if it is from pmc or pubmed and the second one is the article id) in tsv format)")
args = parser.parse_args()
r=Running_config()
r.run(args)
