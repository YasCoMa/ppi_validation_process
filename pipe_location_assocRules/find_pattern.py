from pygosemsim import graph, download
from rdflib import Graph
import urllib.request
import pandas as pd
import os

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori
from mlxtend.frequent_patterns import association_rules

#from apyori import apriori

class Filter_pattern_localization:

    def convert_owl_to_rdf(self):
        link = "http://owl.cs.manchester.ac.uk/converter/convert?ontology=http://purl.obolibrary.org/obo/go.owl&format=RDF/XML"
        f = urllib.request.urlopen(link)
        file = f.read()
        f=open("go.rdf","w")
        f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'').replace('>"','>'))
        f.close()

    def load_go(self):
        got={}
        ide=""
        f=open("go.obo","r")
        for line in f:
            l=line.replace("\n","")
            if(l.find("[Term]")!=-1):
                if(ide!=""):
                    if(class_=="cellular_component"):
                        got[ide]=names
                ide=""
                class_=""
                names=[]
            
            if(l.find("id: ")!=-1 and l.find("_id: ")==-1):
                ide=l.split(": ")[1]

            if(l.find("namespace: ")!=-1):
                class_=l.split(": ")[1]

            if(l.find("name: ")!=-1):
                names.append(l.split(": ")[1].lower())

            if(l.find("synonym: ")!=-1):
                names.append(l.split('"')[1].lower())

            if(l.find("[Typedef]")!=-1):
                if(class_=="cellular_component"):
                    got[ide]=names
                break

        f.close()

        return got

    def find_go_subcellular(self):
        f=open("mapping_uniprotloc_go.tsv","w")
        f.close()
        
        got=self.load_go()

        """graph_go = graph.from_resource("go")

        dict_={}
        for go_term in graph_go.node:
            namespace = graph_go.node[go_term]['namespace']
            if(namespace=="cellular_component"):
                dict_[go_term]=graph_go.node[go_term]['name']"""
        
        f=open("locations_uniprot.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            loc=l[3].lower()
            loc=loc.replace("peroxisome","peroxisomal")
            for go in got:
                
                if(loc in got[go] ):
                    with open("mapping_uniprotloc_go.tsv","a") as gf:
                        gf.write("%s\t%s\t%s\t%s\n" %(l[0], go, loc, "|".join(got[go]) ) )
                    break
                else:
                    name=""
                    for n in got[go]:
                        if(n.find(loc)!=-1):
                            name=n
                            break
                    if(name!=""):
                        with open("mapping_uniprotloc_go.tsv","a") as gf:
                            gf.write("%s\t%s\t%s\t%s\n" %(l[0], go, loc, "|".join(got[go]) ) )
                        break
        f.close()

    def get_rdf(self):
        os.system("rm -r rdf_data")
        os.system("mkdir rdf_data")

        proteins = []
        c=0
        f=open("hint_hs.txt","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split("\t")
                if(not l[0] in proteins):
                    proteins.append(l[0])
                if(not l[1] in proteins):
                    proteins.append(l[1])
            c+=1
        f.close()

        for p in proteins:
            try:
                link = "https://www.uniprot.org/uniprot/"+p+".rdf"
                f = urllib.request.urlopen(link)
                file = f.read()
                f=open("rdf_data/"+p+".rdf","w")
                f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'').replace('>"','>'))
                f.close()
            except:
                pass

    def get_go_annot_from_rdf(self):
        f=open("annotation_data.tsv", "w")
        f.close()
        
        try:
            graph_go = graph.from_resource("go")
        except:
            download.obo("go")
            graph_go = graph.from_resource("go")

        proteins = []
        c=0
        f=open("hint_hs.txt","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split("\t")
                if(not l[0] in proteins):
                    proteins.append(l[0])
                if(not l[1] in proteins):
                    proteins.append(l[1])

            c+=1

        f.close()

        for p in proteins:
            protein=p

            features=[]
            go_cc_ids=[]
            go_mf_ids=[]
            go_bp_ids=[]
            ko_ids=""
            pfam_ids=[]

            try:
                g=Graph()
                g.parse("rdf_data/"+protein+".rdf", format="xml")
                
                results=g.query("""
                SELECT distinct ?o ?c
                WHERE{
                        <http://purl.uniprot.org/uniprot/"""+protein+"""> <http://purl.uniprot.org/core/classifiedWith> ?o .
                        <http://purl.uniprot.org/uniprot/"""+protein+"""> <http://www.w3.org/2000/01/rdf-schema#seeAlso> ?c .
                }
                """)
                namespace=""
                for row in results.result:
                    info=row[0]

                    if(info.find("obo/GO_")!=-1):
                        go=info.replace("http://purl.obolibrary.org/obo/","")
                        go_term=go.replace("_", ":")
                        #namespace = str(self.get_namespace_from_go(go_term.replace(":","_")))
                        if(go_term in graph_go.node):
                            namespace = graph_go.node[go_term]['namespace']

                            if(namespace=="cellular_component"):
                                if(not go_term in go_cc_ids):
                                    go_cc_ids.append(go_term)
                            if(namespace=="molecular_function"):
                                if(not go_term in go_mf_ids):
                                    go_mf_ids.append(go_term)
                            if(namespace=="biological_process"):
                                if(not go_term in go_bp_ids):
                                    go_bp_ids.append(go_term)    
                        #go_ids.append(str(go_term[0]))

                    info2=row[1]
                    if(info2.find("pfam")!=-1 and not(info2.find("/supfam")!=-1)):
                        pfam=info2
                        pfam=pfam.replace("http://purl.uniprot.org/pfam/","")
                        if(not pfam in pfam_ids):
                            pfam_ids.append(pfam)

                    if(info2.find("ko")!=-1):
                        ko=info2
                        ko_ids=ko.replace("http://purl.uniprot.org/ko/","")

                if(len(go_cc_ids)==0):
                    go_cc_ids.append("None")
                if(len(go_mf_ids)==0):
                    go_mf_ids.append("None")
                if(len(go_bp_ids)==0):
                    go_bp_ids.append("None")
                if(len(pfam_ids)==0):
                    pfam_ids.append("None")
                if(ko_ids==""):
                    ko_ids="None"

            except:
                go_cc_ids.append("None")
                go_mf_ids.append("None")
                go_bp_ids.append("None")
                pfam_ids.append("None")
                ko_ids="None" 

            features.append(go_cc_ids)
            features.append(go_mf_ids)
            features.append(go_bp_ids)
            features.append(ko_ids)
            features.append(pfam_ids)

            txt=protein+"\t"
            b=0
            for feature in features:
                if(not(str(type(feature))=="<class 'str'>")):
                    a=0
                    for fea in feature:
                        txt+=str(fea)
                        if(a!=len(feature)-1):
                            txt+=" "
                        a+=1

                    if(b!=len(features)-1):
                        txt+="\t"
                    b+=1
                else:
                    txt+=feature+"\t"
            
            with open("annotation_data.tsv", "a") as g:
                g.write(txt+"\n")

    def load_data_for_itemset_mining(self):
        wl=[]
        f=open("mapping_uniprotloc_go.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            wl.append(l[1])
        f.close()
        #print(wl)

        infop={}
        f=open("annotation_data.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[1]!="" and l[1]!="None"): # 1 for cellular component
                infop[l[0]]=l[1].split(" ")
            else:
                infop[l[0]]=[]

        f.close()

        dataset=[]
        c=0
        f=open("hint_hs.txt","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split("\t")
                #if(l[0] in wl and l[1] in wl):
                if(len(infop[l[0]])>0 and len(infop[l[1]])>0 ):
                    ds=[]
                    for a in infop[l[0]]:
                        ds.append(a)
                    for a in infop[l[1]]:
                        ds.append(a)
                    
                    dataset.append(ds)
            c+=1
        f.close()

        return dataset

    def evaluate(self, threshold):
        
        wl={}
        nl={}
        f=open("mapping_uniprotloc_go.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            wl[l[0]]=l[1]
            nl[l[0]]=l[2]
        f.close()

        c=0
        missing=[]
        f=open("missing_locations.csv","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split("\t")
                if(l[0] in wl.keys()):
                    missing.append(wl[l[0]])
                else:
                    print(l[0],"\t",nl[l[0]])
            c+=1
        f.close()

        dic={}
        for m in missing:
            dic[m]=0

        c=0
        f=open("rules.csv","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").replace('"',"").replace('frozenset(',"").replace(')',"").replace('{',"").replace('}',"").replace("'","").replace(' ',"").split("\t")
                all=[]
                ant = l[0].split(",")
                for a in ant:
                    if( not a in all):
                        all.append(a)

                con = l[1].split(",")
                for a in con:
                    if( not a in all):
                        all.append(a)

                for el in all:
                    if(el in dic.keys()):
                        dic[el]+=1
            c+=1
        f.close()

        count_coverage=0
        for el in dic.keys():
            if(dic[el]>0):
                print(el)
                count_coverage+=1

        coverage=count_coverage/len(missing)
        print(count_coverage, len(missing))

        return coverage>=threshold


    def itemset_pattern(self):
        graph_go = graph.from_resource("go")

        dataset=self.load_data_for_itemset_mining()
       
        """dataset = [['Milk', 'Onion', 'Nutmeg', 'Kidney Beans', 'Eggs', 'Yogurt'],
           ['Dill', 'Onion', 'Nutmeg', 'Kidney Beans', 'Eggs', 'Yogurt'],
           ['Milk', 'Apple', 'Kidney Beans', 'Eggs'],
           ['Milk', 'Unicorn', 'Corn', 'Kidney Beans', 'Yogurt'],
           ['Corn', 'Onion', 'Onion', 'Kidney Beans', 'Ice cream', 'Eggs']]"""

        te = TransactionEncoder()
        te_ary = te.fit(dataset).transform(dataset)
        df = pd.DataFrame(te_ary, columns=te.columns_)

        flag_reach_coverage=False
        i=0.3
        while(i>=0.01 and not flag_reach_coverage):
            frequent_itemsets = apriori(df, min_support=i, use_colnames=True)
            #rules = association_rules(frequent_itemsets, metric="lift", min_threshold=1.2)
            rules = association_rules(frequent_itemsets, metric="confidence", min_threshold=i)
            print("--->",i, len(rules))
            rules.to_csv(r'rules.csv', sep="\t", index = False)
            flag_reach_coverage=self.evaluate(0.6)
            #print(i, flag_reach_coverage)
            i-=0.01

        f=open("rules_go.tsv","w")
        f.write("antecedente\tconsequente\tsupport\tconfidence\tlift\n")
        f.close()
        
        f=open("rules_name.tsv","w")
        f.write("antecedente\tconsequente\tsupport\tconfidence\tlift\n")
        f.close()

        got=self.load_go()

        c=0
        f=open("rules.csv","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").replace('"',"").replace('frozenset(',"").replace(')',"").replace('{',"").replace('}',"").replace("'","").replace(' ',"").split("\t")
                
                ant = l[0].split(",")
                ant_name=[]
                for a in ant:
                    ant_name.append(got[a][0])

                con = l[1].split(",")
                con_name=[]
                for a in con:
                    con_name.append(got[a][0])

                with open("rules_name.tsv","a") as gf:
                    gf.write("%s\t%s\t%s\t%s\t%s\n" %( ";".join(ant_name), ";".join(con_name), l[4], l[5], l[6] ) )

                with open("rules_go.tsv","a") as gf:
                    gf.write("%s\t%s\t%s\t%s\t%s\n" %( ";".join(ant), ";".join(con), l[4], l[5], l[6] ) )
            
            c+=1
        f.close()

        """
        f=open("rules_go.tsv","w")
        f.write("antecedente\tconsequente\tsupport\tconfidence\tlift\n")
        f.close()
        
        f=open("rules_name.tsv","w")
        f.write("antecedente\tconsequente\tsupport\tconfidence\tlift\n")
        f.close()
        c=0
        for r in list(rules):
            if(c>0):
                print(r)
                ant = list(r[1])
                print(ant)
                ant_name=[]
                for a in ant:
                    ant_name.append(graph_go.node[a]['name'])

                con = list(r[2])
                print(con)
                con_name=[]
                for a in ant:
                    con_name.append(graph_go.node[a]['name'])
                
                with open("rules_name.tsv","a") as gf:
                    gf.write("%s\t%s\t%s\t%s\t%s\n" %( ";".join(ant_name), ";".join(con_name), r[5], r[6], r[7] ) )

                with open("rules_go.tsv","a") as gf:
                    gf.write("%s\t%s\t%s\t%s\t%s\n" %( ";".join(ant), ";".join(con), r[5], r[6], r[7] ) )
            c+=1"""
        #rules.to_csv(r'rules.csv', index = False)

        """association_rules = apriori(dataset, min_support=0.0045, min_confidence=0.2, min_lift=3, min_length=2)
        association_results = list(association_rules)
        #association_rules.to_csv(r'rules.csv', index = False)
        for item in association_rules:
            # first index of the inner list
            # Contains base item and add item
            pair = item[0] 
            items = [x for x in pair]
            print("Rule: " + items[0] + " -> " + items[1])

            #second index of the inner list
            print("Support: " + str(item[1]))

            #third index of the list located at 0th
            #of the third index of the inner list

            print("Confidence: " + str(item[2][0][2]))
            print("Lift: " + str(item[2][0][3]))
            print("=====================================") """

# http://rasbt.github.io/mlxtend/user_guide/frequent_patterns/association_rules/#example-2-rule-generation-and-selection-criteria

class Filter_by_association_rules:
    def mapping(self):
        mapp={}
        f=open("mapping_hgnc_uniprot.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            mapp[l[0]]=l[1]
        f.close()

        return mapp

    def mapping_reverse(self):
        mapp={}
        f=open("mapping_hgnc_uniprot.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            mapp[l[1]]=l[0]
        f.close()

        return mapp

    def load_rules(self):
        rules=[]
        f=open("rules_go.tsv")
        n=0
        for line in f:
            if(n>0):
                l=line.replace("\n","").split("\t")
                rules.append([l[0].split(";"), l[1].split(";")])
            n+=1
        f.close()

        return rules
    
    def load_annotations(self, proteins):
        infop={}
        for p in proteins:
            if(os.path.isfile("annotation_data/"+p+".tsv")):
                f=open("annotation_data/"+p+".tsv","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    if(l[1]!="None" and l[1]!=""):
                        ccs=l[1].split(" ")
                        infop[l[0]]=ccs
                    else:
                       infop[l[0]]=[]
                f.close()
            else:
                infop[p]=[]

        return infop

    def run(self, folder, file_pairs):
        print("Loading rules...")
        rules=self.load_rules()

        print("Loading pairs...")
        pairs=[]
        f=open(folder+file_pairs,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            pr=[l[0], l[1] ]
            if(not pr in pairs):
                pairs.append(pr)
        f.close()

        proteins=[]
        for p in pairs:
            if(not p[0] in proteins):
                proteins.append(p[0])
            if(not p[1] in proteins):
                proteins.append(p[1])

        print("Loading annotation about locations....")
        infop=self.load_annotations(proteins)

        print("Filtering pairs according to the rules...")
        nrules=[]
        filtered_proteins=[]
        n=0
        for p in pairs:
            c1=0
            c2=0

            nr=0
            for r in rules:
                if(n==0):
                    nrules.append([])
                
                if( len(infop[p[0]]) > 0 ):
                    for cc in infop[p[0]]:
                        if(cc in r[0]):
                            c1+=1

                if( len(infop[p[1]]) > 0 ):
                    for cc in infop[p[1]]:
                        if(cc in r[1]):
                            c2+=1
                    
                if(c1>0 and c2>0):
                    nrules[nr].append(p)
                nr+=1
            
            if(c1>0 and c2>0):
                filtered_proteins.append(p)
            n+=1
        
        print("Exporting filtered pairs...")
        f=open(folder+"filtered_pairs_by_assocrules.tsv","w")
        f.close()
        mapp=self.mapping_reverse()
        for p in filtered_proteins:
            with open(folder+"filtered_pairs_by_assocrules.tsv","a") as fg:
                fg.write("%s|%s\t%s|%s\t1\n" %( p[0], mapp[p[0]], p[1], mapp[p[1]] ) )

        print("File containing the filtered pairs can be found in filtered_pairs_by_assocrules.tsv")

        print("Exporting report associated to the rules...")
        f=open(folder+"report_rules.tsv","w")
        f.write("antecedente\tconsequente\tnumero de pares\n")
        f.close()
        # reporting rules
        nr=0
        for r in rules:
            with open(folder+"report_rules.tsv","a") as fg:
                fg.write("%s\t%s\t%i\n" %(";".join(r[0]),";".join(r[1]), len(nrules[nr])) )
            nr+=1

        print("Report about rules participation can be found in report_rules.tsv")

# Other functions to discover the rules:
#a=Filter_pattern_localization()
#a.get_rdf()
#a.get_go_annot_from_rdf()
#a.itemset_pattern()

class Running_config:
    def run(self, args):
        if(args.folder!="" and args.interactome_file!="" ):
            if(os.path.isdir(args.folder) and os.path.isfile(args.folder+args.interactome_file) ):
                a=Filter_by_association_rules()
                a.run(args.folder, args.interactome_file)
            else:
                print("Error: There are invalid folder or files, some of them were not found")
        else:
            print("Error: There are missing parameters")

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Pipeline to filter pairs according to association rules extracted from HINT database', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-fo", "--folder", action="store", help="Required - Folder to store the files (use the folder where the other required file can be found)")
    parser.add_argument("-if", "--interactome_file", action="store", help="Required - File with the pairs (two columns with uniprot identifiers in tsv format)")
    args = parser.parse_args()
    r=Running_config()
    r.run(args)
