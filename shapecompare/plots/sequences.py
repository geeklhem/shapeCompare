palette = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
model_palette = ["#FF0000",'#66CCFF'] 
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter

def plot(plt,shapeData,display="all",models={}):	
    """construct the sequence plot"""
    if display == "all":
          display = [1] * len(shapeData)

    patches = []
    n = -1
    line = 0
    for k, sd in shapeData.iteritems(): #For each type of RNA

        nb_exp =  sd.reactivity.shape[0]
        
        #Find if there is any model related to this RNA type and list them
        selected_models = []
        for j, m in models.iteritems():
            if m.rna == k:
                selected_models.append(m)
        nb_models = len(selected_models)
        
        for (i,base) in enumerate(sd.sequence_alg):
            if base != "-": #Alignement
                n += 1

                # Experimental data
                for exp in range(sd.reactivity.shape[0]):
                    plt.text(i+0.5,nb_models+exp+line+0.4,str(base))
                    if n%10==0:
                        plt.text(i,nb_models+exp+line+0.6,str(n))
                    color = int(sd.reactivities_class[exp,n])
                    patches.append(mpatches.Rectangle((i,nb_models+exp+line),
                                                      1,
                                                      1,
                                                      edgecolor='none',
                                                      fc=palette[color]))
       
                # Display models associated with this typr of RNA (if any)
                for j , m in enumerate(selected_models):
                    color = m.app[n]
                    plt.text(i+0.5,j+line+0.4,str(base))
                    if n%10==0:
                        plt.text(i,j+line+0.6,str(n))
                  
                    patches.append(mpatches.Rectangle((i,j+line),
                                                      1,
                                                      1,
                                                      edgecolor='none',
                                                      fc=model_palette[color]))
                    
        line += nb_models + nb_exp

        #Label lines
        for w, name in enumerate(sd.names):
                plt.text(i-1,line-nb_exp+w+0.1,"["+k+"]: "+name)
        for j , m in enumerate(selected_models):
                plt.text(i-1,line-1-w-nb_models+j+0.1,"["+k+"]: "+m.name+" model")
        n = -1

    collection = PatchCollection(patches,match_original=True)
    plt.add_collection(collection)
    plt.axis([i,0,0,line])
    plt.xaxis.set_major_formatter(NullFormatter())
    plt.yaxis.set_major_formatter(NullFormatter())
