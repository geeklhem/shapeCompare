from numpy import *
import options
from fasta import Fasta

class Model():
    """Encapsulate a secondary structure model informations."""
    sequences = Fasta(options.FASTA_FILE)
    def __init__(self,name):
        self.name = name

        try:
            self.rna= str(recfromtxt(options.MODEL_PATH+self.name+"_seq.txt"))
        except:
            self.rna = self.sequences.sequences.keys()[-1] #by default sequence is the first of the fasta file


        try:
            self.position = recfromtxt(options.MODEL_PATH+
                                       self.name+"_pos.csv",
                                       delimiter=",",dtype=float)
       
            self.position = self.position/amax(self.position)
      
        except:
            print("Position file of "+name+" model absent or corrupted")
            self.position = zeros((300,2))
       
        try:
            self.app = recfromtxt(options.MODEL_PATH+self.name+"_app.csv")
        except:
            print("Appariment file of "+name+" model absent or corrupted")
            self.app = zeros(300)

        def __str__(self):
            '''String representaion of objects'''
            return(self.name+" model")
  
