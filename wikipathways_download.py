import pywikipathways as pwpw
import os


hs_pathways = pwpw.list_pathways('Homo sapiens')

if not os.path.exists('KGML'):
    os.mkdir('KGML')

for id in hs_pathways.id.tolist():
    gpml = pwpw.get_pathway(id)
    if gpml is not None:
        f=open("KGML/"+id+".gpml","w")
        f.write(gpml)
        f.close()