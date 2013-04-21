from numpy import *
import scipy.stats as st
from fasta import Fasta
from random import random
import options

def normValue(shape_data,column):
    av = 0
    shape_data = shape_data[shape_data[:,2].argsort()]
    lg = shape(shape_data)[0]
    rg = range(int(ceil((lg/100.)*2)),int(ceil((lg/100.)*10)))
    for i in rg:
        av += shape_data[-i,2]
    av = av / len(rg)
    return(av)


class ShapeData():
    '''Contain peak data extracted from multiple shape convertor object, align them to the sequence file given and perform statistical operations'''
    alignement = Fasta(options.FASTA_ALG)
    sequences = Fasta(options.FASTA_FILE)

    def __init__(self,data,sequenceFile):
        self.sequence = self.sequences.sequences[sequenceFile]
        self.sequence_alg = self.alignement.sequences[sequenceFile]
        
        self.names= []
        self.nb_files = len(data)

        # ALIGN DATA ON SEQUENCE
        self.reactivity = zeros((self.nb_files,len(self.sequence)))
        self.reactivities_class = zeros((self.nb_files,len(self.sequence)))
        self.reactivities_class += 5
        for k,d in enumerate(data):
            len_experiment = len(d.peak_reactivity)
            for j in range(len_experiment):
                pos = -(j + d.offset_bases)
                self.reactivity[k,pos] = d.peak_reactivity[j]
                if d.colors[j] == -1:
                    self.reactivities_class[k,pos] = 5
                else:
                    self.reactivities_class[k,pos] = d.colors[j]
            self.names.append(d.name)
        
        ## STATISTICAL OPERATION
        self.mean_flex = mean(self.reactivity,0)
        self.var = std(self.reactivity,0)
        self.mean_reactivity_class = zeros(len(self.sequence))
        self.variability_class = zeros(len(self.var))


        if self.nb_files>1:
            self.compute_variability()

    def compute_cor(self,base):
        r = []
        for k in range(len(self.reactivity[0,:])):
            if self.variability_class[k] >= 1:
                coef = st.pearsonr(x=self.reactivity[:,base],y=self.reactivity[:,k])
                r.append([coef[0],coef[1],k])
        r = array(r)
        r= r[r[:,1].argsort()]
        return(r)
        

    def __str__(self):
        '''String representaion of objects'''
        return("<Shape Variability Object on "+str(self.nb_files)+" file(s)>")
    
    def compute_variability(self):
       

        for i, flex in enumerate(self.mean_flex):
            if flex<0: #Unknown
                self.mean_reactivity_class[i] = 0
            elif flex>=1: #Unknown -> Grey
                self.mean_reactivity_class[i] = 5
            else:
                self.mean_reactivity_class[i] = int(floor(flex*5)) 
        
        #I use var to name variable but i've ended in using standard deviation instead of variance.
        self.varmean = mean(self.var)
        self.sd_of_var = std(self.var)

        for i,v in enumerate(self.var):
            self.variability_class[i] = int(abs(self.varmean - v)/self.sd_of_var)
            if self.variability_class[i] >= 1:
                self.mean_reactivity_class[i] = 6
                

class ShapeConvertor():
    '''Convert data from *.shape files created by shapeFinder'''
    sequences = Fasta(options.FASTA_FILE)
  

    def __init__(self,path,name):

        # LOADING
        self.name = name
        self.path = path + '/'
        self.load(self.path)
              
        # SEQUENCE RECOGNITION
        try:
            self.seq_type = str(recfromtxt(self.path+'seq.txt'))
        except:
            self.seq_type = self.sequences.sequences.keys()[-1] #by default sequence is the first of the fasta file
        self.seq = self.sequences.sequences[self.seq_type]        
        self.seq = self.seq[::-1] #Because the data have a reversed sequence.

        # NORMALISE DATA
        self.normalise_seq((2,3))
        self.normalise_shape((0,1))
    
        # SEQUENCE ALIGNEMENT
        self.match_try = 0
        self.match = 0.0
        
        while (self.match<=options.MATCH_THRESHOLD): 
            self.offset = 0
            self.offset_bases = 0
            self.match_try += 1

            # Create an anchor at a random position in the sequence
            size =  10 - int(self.match_try/2)
            self.a_pos = int(random()*(len(self.seq)-size))
            self.a = ""

            for i in range(size):
                    self.a += self.seq[self.a_pos+i]

            #If it can be found, compute the match score.
            if self.compute_offset(self.a,self.a_pos):
                for i,l in enumerate(self.algBases):
                        if l[1] == self.seq[i+self.offset_bases]:
                            self.match +=  1 
                self.match /= len(self.algBases)
        self.match = int(self.match*100) 
        
    def __str__(self):
        '''String representaion of objects'''
        return("<Shape Data Object of \""+self.path+"\" file>")
    
    def compute_offset(self, anchor, anchor_pos):
        '''Compute the offset using the Rabin-Karp algorithm'''
        h_anchor = hash(anchor)
        anchor_l = len(anchor)
        test_len = len(self.seqBases) - anchor_l
        window = []
        match = False
        for i in range(anchor_l):
            window.append(self.algBases[i][1])
        for j in range(i+1,test_len):
            window.pop(0)
            window.append(self.algBases[j][1])
            if h_anchor == hash(''.join(window)) :
                self.offset = 2000 - self.algBases[j+1-anchor_l][0]
                self.offset_bases = anchor_pos-(j+1-anchor_l)
                match = True
                break
        return match
    

    def load(self,path):
        '''load data from the .shape file'''
        #Trace : size = (NbOfpx,4), value for each channel 
        self.rawTrace = recfromtxt(path+'rawTrace')
        self.procTrace = recfromtxt(path+'procTrace')
        self.x = arange(shape(self.procTrace)[0])
        #Sequenced Bases : size = (seqB,2), position and nucleotid of each sequenced bases
        self.seqBases  = recfromtxt(path+'procTrace.aln', skip_header=1,usecols=(0,1))
        self.seqBases.dtype.names = ('px','base')
         
        #Peak : size = (nbOfPeaks,2), value = position (in px) of each peak
        self.peak_pos = recfromtxt(path+'procTrace.pkl')
        self.nb_peaks,self.channels=shape(self.peak_pos)
        self.peak_values = zeros((self.nb_peaks,self.channels))
        self.colors = zeros(self.nb_peaks)
        self.peak_reactivity = zeros(self.nb_peaks)
       
        for c in range(self.channels):
            for p in range(self.nb_peaks):
                self.peak_values[p,c] = self.procTrace[self.peak_pos[p,c],c]
        
        #Aligned Bases size = (nbOfPeaks,2), value = position in px and nucleotid of aligned bases
        self.algBases = recfromtxt(path+'procTrace.bsl', skip_header=1,usecols=(0,1))
        self.algBases.dtype.names = ('px','base')
       

    def normalise_seq(self,chans):
        '''Normalise data 
        Set the highest peak to 1
        '''
        norm = max(normValue(self.peak_values,chans[0]),normValue(self.peak_values,chans[0]))
        for i in range(shape(self.procTrace)[0]):
            self.procTrace[i,chans[0]] /= norm              
            self.procTrace[i,chans[1]] /= norm
    
    def normalise_shape(self,chans):
        '''1) Exclude the 10% shape- highest
           2) Compute reactivity (shape+)-(shape-)
           3) Set the average reactivity of the 10% highest excluding the first 2% to 1. 
           4) Colorise by reactivity.
        '''
        one = normValue(self.peak_values,chans[0])

        #Shape data : pos shape+, position shap-, valeur shape +, valeur shape -
        shape_data = hstack((self.peak_pos[:,chans], self.peak_values[:,chans]))
        #sort by shape- and discard the 10% strongest
        shape_data = shape_data[shape_data[:,3].argsort()]
        for i in range(int(ceil((self.nb_peaks/100.)*10))):
            shape_data[-i,3] = 10**9
            
        #compute reactivity 
        for i in range(self.nb_peaks):
            shape_data[i,2] = shape_data[i,2] - shape_data[i,3]

        #sort by reactivity and compute average on the 10% strongest
        #except the 2% strongest. 
        av = 0
        shape_data = shape_data[shape_data[:,2].argsort()]
        for i in range(int(ceil((self.nb_peaks/100.)*2)),int(ceil((self.nb_peaks/100.)*10))):
            av += shape_data[-i,2]
        av = av / i
            
        #Colorize
        for i in range(self.nb_peaks):
            
            #Compute peak reactivity
            self.peak_reactivity[i] = self.procTrace[self.peak_pos[i,chans[0]],chans[0]] 
            self.peak_reactivity[i] -= self.procTrace[self.peak_pos[i,chans[1]],chans[1]]
            self.peak_reactivity[i] /= av           
       
            #Assign category
            if self.peak_reactivity[i]<0:
                self.colors[i] = 0
            elif self.peak_reactivity[i]>=1:
                self.colors[i] = 4
            else:
                self.colors[i] = int(self.peak_reactivity[i]*5) 
            
            if self.procTrace[self.peak_pos[i,chans[1]],chans[1]]/one>0.90:
                    self.colors[i] = -1

       #Normalise
        for i in range(shape(self.procTrace)[0]):
            self.procTrace[i,chans[0]] /= one             
            self.procTrace[i,chans[1]] /= one
            if self.procTrace[i,chans[1]] > 0.90:
                self.procTrace[i,chans[1]] = 10**9/one
 



if __name__ == "__main__":
    d = ShapeConvertor("data/7SKWT_1m7v_MC6_alin.shape")
    
