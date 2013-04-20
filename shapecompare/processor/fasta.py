import numpy as np

class Fasta():
    def __init__(self,path):
        data = np.recfromtxt(path)
        i = -1
        self.sequences = {}
        #self.names={}
        for l in data:
            if ">" in l:
                name = l[1:] #To strip the >
                self.sequences[name]=('')
            else:
                self.sequences[name]+=l

if __name__ == "__main__":
    a = Fasta("../../data/sequences/7SK_alg.fasta")
    print(a.sequences)
