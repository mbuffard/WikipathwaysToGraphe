from xml.dom import minidom
import converter
import os
from os import listdir
compteur_noeuds_none=0
compteur_noeuds_ok=0
compteur_pathway_vide=0
compteur_pathway=0






for pathway_file in listdir('GPML/'):
    
    pathway_name=pathway_file.split('.')[0]
    f=minidom.parse('GPML/'+pathway_file)    
    node_list={}
    #avoid nodes that doesn't correspond to proteins, Type de l'instance DataNode
    NODES_IGNORED=['Undefined','DNA','RNA','Complex','Metabolite','Pathway','Disease','Phenotype','Alias','Event']
    DATABASE_ok=['Uniprot-TrEMBL','Uniprot-SwissProt']


    #look at every node and translate name into unique uniprot identifier
    node_file=open("pathways_nodes/"+pathway_name+"_nodes.txt","w")
    for datanode in f.getElementsByTagName('DataNode'):
        if 'Type' not in datanode.attributes or 'GraphId' not in datanode.attributes:
            continue
        if 'GroupRef' not in datanode.attributes and datanode.attributes['Type'].value not in NODES_IGNORED:
            name = converter.handler.import_symbol(datanode.attributes['TextLabel'].value)
            #map unrecognized node by id and database
            if name is None:  
                for xref in datanode.getElementsByTagName('Xref'):
                    db_name=xref.attributes['Database'].value
                    if db_name in DATABASE_ok:
                        id=xref.attributes['ID'].value
                        clean_id=converter.handler.clean_uid(id)
                        
                        if clean_id is not None  :
                            clean_name=converter.handler.to_symbol(clean_id)
                            
                            node_list.update({datanode.attributes['GraphId'].value:clean_id})
                            
                            node_file.write(clean_id+"\t"+str(clean_name)+"\n")
                            compteur_noeuds_ok+=1
                        else:
                            compteur_noeuds_none+=1

                        #print("database ok "+converter.handler.to_symbol(clean_id))
                    #print(converter.handler.import_symbol('PIK3R5'))
                    
                
            else:
                compteur_noeuds_ok+=1
                if 'GraphId' not in datanode.attributes:
                    continue
                node_list.update({datanode.attributes['GraphId'].value:name})
                node_file.write(name+"\t"+datanode.attributes['TextLabel'].value+"\n")


    #extract interactions and write a source target tabulated file
    if not os.path.exists('pathways'):
        os.mkdir('pathways')

    g=open('pathways/'+pathway_name+'.txt',"w")
    i=0
    interaction_dico={}

    for interaction in f.getElementsByTagName('Interaction'):
        source=''
        target=''
        i+=1
        for Graph in interaction.getElementsByTagName('Graphics'):
            interact_group=0
            
            for point in Graph.getElementsByTagName('Point'):
                #if the node has been translated into uniprot id
                if 'GraphRef' in point.attributes:
                    if point.attributes['GraphRef'].value in node_list:
                        #node with ArrowHead are the targets
                        if 'ArrowHead' in point.attributes:
                            target = node_list[point.attributes['GraphRef'].value]
                        else:
                            source = node_list[point.attributes['GraphRef'].value]
                        #print (source+"->"+target)
                #the  node is a group or the HGNC is not good
            if source!='' and target !='' and source!=target :
            
                if source in interaction_dico:
                    if target not in interaction_dico[source]:
                        interaction_dico[source].append(target)
                else:
                    interaction_dico[source]=[target]
    for key in interaction_dico:
        for value in interaction_dico[key]:
            g.write(key+"\t"+value+"\t\t\n")
    g.close()
    if interaction_dico:
        compteur_pathway+=1
    else:
        compteur_pathway_vide+=1
        
#print (interaction_dico)
#print(node_list)
print('nombre de noeuds non traduits')
print((compteur_noeuds_none))
print('nombre de noeuds traduits')
print((compteur_noeuds_ok))
print('nombre de pathways transformés :')
print (compteur_pathway)
print('nombre de pathways non transformés :')
print(compteur_pathway_vide)

#Supress all the empty pathway files
for file in listdir('pathways/'):
    if os.path.getsize('pathways/'+file) == 0:
        os.remove('pathways/'+file)
        node_file_name=file.split('.')[0]+"_nodes.txt"
        os.remove('pathways_nodes/'+node_file_name)



# Write the member file used by enriched analysis
member_file=open("pathways/members.txt","w")
for pathway_file in listdir('pathways/'):
    if pathway_file=="members.txt":continue
    f=minidom.parse('GPML/'+pathway_file.split('.')[0]+".gpml")
    
    print(f.attributes)
    pathway_name=pathway_file.split('.')[0]
    member_file.write(pathway_name+"\t")
    if not f.attributes:
        member_file.write(pathway_name+"\t")
    else:
        if 'Pathway' not in f.attributes:
            member_file.write("no description"+"\t")
        else:
            if "Name" in f.attributes['Pathway']:
                name=f.attributes['Pathway'].attributes["Name"].value
                member_file.write(name+"\t")
    node_file=open("pathways_nodes/"+pathway_file.split(".")[0]+"_nodes.txt")

    for line in node_file:
        member_file.write(line.strip().split("\t")[0]+",")
    member_file.write("\n")
member_file.close()


