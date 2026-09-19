from Bio.SeqUtils import ProtParam
from Bio import SeqIO
import requests
from Bio.SeqUtils.ProtParamData import kd,Flex
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+gene:HSPA1A+AND+organism_name:Homo+sapiens&format=fasta"
data = requests.get(url).text
with open ("HSPA1A.fasta","w") as f :
    f.write(data)
parser = SeqIO.parse("HSPA1A.fasta","fasta")
for record in parser :
  sequence = str(record.seq)
  analysis = ProtParam.ProteinAnalysis(sequence)
  data_stored= {}
  data_stored["no_of_aa"] = len(sequence)
  data_stored["freq of each aa"] = analysis.count_amino_acids()
  data_stored["weight"]= analysis.molecular_weight()
  data_stored["aromaticity"] = analysis.aromaticity()
  data_stored["instability_index"] = analysis.instability_index()# below 40- stable , above 40 - unstable
  data_stored["sec_str_fraction"]= analysis.secondary_structure_fraction()# tuple(helix,turn,sheet)
  data_stored["hydrophabicity index"]= analysis.gravy() # positive value- hydrophobic , negative value - hydrophilic
  hydrophobicity = np.array(analysis.protein_scale(kd, 9)) # kd - hydrophobicity of 9 sliding windows
  flexibility = np.array(analysis.protein_scale(Flex, 9)) # flex - flexibility of 9 sliding windows
  data_stored["regions"] = {
    "NBD": {                                     # NBD- nucleotide binding region
        "start": 2,
        "end": 386,
        "description": "Nucleotide-binding domain"
    },

    "SBD": {                                    # SBD - substrate binding region
        "start": 394,
        "end": 509,
        "description": "Substrate-binding domain"
    },

    "glycine_rich_region": {
        "start": 614,
        "end": 633,
        "description": "Glycine-rich region"
    },

    "disordered_region": {
        "start": 614,
        "end": 641,
        "description": "Disordered region"
    }
  }
  df = pd.DataFrame(
     { "position" : np.arange(5,638),
        "hydrophobicity" : hydrophobicity,
        "flexibility" : flexibility }
  )
  df.to_csv("HSPA1A_data.csv", index=False)
  data_stored["corelation between hydrophobicity and flexibility"] = df.corr()
  plt.plot(df["position"],df["hydrophobicity"],color= "green", label= "hydrophobicity")
  plt.plot(df["position"],df["flexibility"],color= "red", label= "flexibility")
  plt.xlabel('position')
  plt.ylabel('hydrophobicity and flexibility corelation')
  plt.title('Project Comparison')
  plt.grid(True) # Adds a background grid for better metric evaluation
  plt.legend()   # Displays the labels assigned in the plot function
  
  plt.show()

  print(data_stored)
