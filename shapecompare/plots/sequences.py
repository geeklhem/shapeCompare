palette = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter

def plot(plt,shapeData):	
		"""construct the sequence plot"""
		patches = []
		react_color = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
		n = 0
		exp_offset = 0
		for k, sd in shapeData.iteritems():
			for (i,base) in enumerate(sd.sequence_alg):
				if base != "-":
					for exp in range(sd.reactivity.shape[0]):
						plt.text(i+0.5,exp+exp_offset+0.4,str(base))
						if n%10==0:
							plt.text(i,exp+exp_offset+0.6,str(n))
						color = int(sd.reactivities_class[exp,n])
						patches.append(mpatches.Rectangle((i,exp+exp_offset),
										  1,
										  1,
										  edgecolor='none',
										  fc=react_color[color]))
					n = n+1
		                
			exp_offset += sd.reactivity.shape[0]
			n = 0

		collection = PatchCollection(patches,match_original=True)
		plt.add_collection(collection)
		plt.axis([i,0,0,exp_offset])
                plt.xaxis.set_major_formatter(NullFormatter())
                plt.yaxis.set_major_formatter(NullFormatter())
                
		
