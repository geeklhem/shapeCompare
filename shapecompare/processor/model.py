from numpy import *

class Model():
    """Encapsulate a secondary structure model informations."""
    
    def __init__(self,name):
        self.name = name
        try:
            self.position = recfromtxt("../data/models/"+self.name+"_pos.csv",delimiter=",",dtype=float)
       
            self.position = self.position/amax(self.position)
            #self.position = self.position[::-1]
      
        except:
            print("Position file of "+name+" model absent or corrupted")
            self.position = zeros((300,2))
       
        try:
            self.app = recfromtxt("../data/models/"+self.name+"_app.csv")
        except:
            print("Appariment file of "+name+" model absent or corrupted")
            self.app = zeros(300)

       
 #       print(self.position)
 #       print(self.app)

        def __str__(self):
            '''String representaion of objects'''
            return(self.name+" model")
  
