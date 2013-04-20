try:
    	from traits.api import *
	from traitsui.api import *	                    
except ImportError:
	from enthought.traits.api import *
	from enthought.traits.ui.api import *
	
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter() #Nolabels
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from numpy import empty
from shapeprocessor.model import Model

class ModelPanel(HasTraits):
	"""Model Panel"""
	dview = Property
	plot_button = Button("Plot")
	display_sequences = Bool(True)
	display_model2d = Bool()
	display_meanReactivity = Bool(True)
	display_sdReactivity = Bool()
	display_sbv = Bool()
	selected_base = Int(237)
	seq_names = Str("Order:"+"\n"*15)
	selected_model = Int(0)
	select_model2d = Str("7SK_Human")
	select_reactivity = Str("7SK_Human")
	select_sdreactivity = Str("7SK_Human")

	def __init__(self,panel,figure,data):
		self.panel = panel
		self.figure = figure
		self.data = data
		self.shapeData = self.panel.shapeData
		self.display = (self.display_sequences,self.display_model2d,self.display_meanReactivity)  
		self.models = []
		self.models.append(Model("wassarman"))
		self.models.append(Model("marz"))
		self.models.append(Model("eilebrecht"))
		self.model_names = {}
		

		for i,m in enumerate(self.models):
			self.model_names[i] = m.name

	def _get_dview(self):
		trait_names = self.editable_traits()
		trait_names.remove('dview')
		return(View(Group(Item("display_sequences",label="Display:"),
				  Item("seq_names", show_label=False,style="readonly"),
				  label = "Sequences Plot",
				  show_border = True),
			    Group(Item("display_model2d",label="Display:"),
				  Item("selected_model", 
				       editor=EnumEditor(values=self.model_names),),
				  Item(name='select_model2d',
				       label="RNA",
				       editor=EnumEditor(values=self.shapeData.keys(),cols=2)), 
				  label = "2D Sequences Plot",
				  show_border = True),
			    Group(Item("display_meanReactivity",label="Display:"),
				  Item(name='select_reactivity',
				       label="RNA",
				       editor=EnumEditor(values=self.shapeData.keys(),cols=2)),
				  label = "Mean reactivity Plot", 
				  show_border = True),
			    Group(Item("display_sdReactivity",label="Display:"),
				  Item(name='select_sdreactivity',
				       label="RNA",
				       editor=EnumEditor(values=self.shapeData.keys(),cols=2)), 
				  label = "Standard deviation of reactivity Plot",
				  show_border = True),
			    Item("plot_button",show_label=False),
			    handler=ModelHandler(self)))

	def data_loaded(self):
		self.seq_names = "Order:\n"
		for k, sd in self.shapeData.iteritems():
			self.seq_names += "\n  " + k + " : \n"
			for n in reversed(sd.names):
				self.seq_names += n + "\n" 
		self.trait_property_changed('dview',None,self.dview)

	def _plot_button_fired(self):
                #Clear window
		self.figure.clf()
		self.plots = []
		self.plot_type = []

		self.display = (self.display_sequences,self.display_model2d,self.display_meanReactivity,self.display_sdReactivity,self.display_sbv)  
		nb_of_plots = 0
		i = 0
		for plot in self.display:
			if plot: 
				nb_of_plots += 1		
		for j, plot in enumerate(self.display):
			if plot: 
				if j == 0:
					self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1))
					self.sequences_plot(self.plots[i])
				if j == 1:
					self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1))
					self.model2d_plot(self.plots[i])
				if j == 2:
					if self.display_sequences:
						self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1,sharex=self.plots[0]))
					else:
						self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1))		
					self.reactivity_plot(self.plots[i])
				if j == 3:
					if self.display_sequences:
						self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1,sharex=self.plots[0]))
					else:
						self.plots.append(self.figure.add_subplot(nb_of_plots,1,i+1))		
					self.sdreactivity_plot(self.plots[i])
				i += 1
			

		for plt in self.plots:
			#plt.axis([0,8000,0,1])
			plt.xaxis.set_major_formatter(nullfmt)
			plt.yaxis.set_major_formatter(nullfmt)
			self.figure.subplots_adjust(hspace=0,bottom=0,top=1,right=1,left=0)
		
                #Draw them
		
		self.figure.canvas.draw()


	def model2d_plot(self,plt):
		react_color = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
		
		app = 0
		for n,(xpos,ypos) in enumerate(self.models[self.selected_model].position):
			print("RNA :"+self.select_model2d)
			color = int(self.shapeData[self.select_model2d].mean_reactivity_class[n])
			if self.models[self.selected_model].app[n]:
				app = "bold"
			else:
				app = "normal"
			plt.text(xpos,ypos,
				 str(self.shapeData[self.select_model2d].sequence[n]),
				 backgroundcolor=react_color[color],
				 weight=app)
			plt.axis([0,1,1,0])
			plt.axis('equal')
		

	def sequences_plot(self,plt):
		"""construct the sequence plot"""
		patches = []
		react_color = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
		n = 0
		exp_offset = 0
		for k, sd in self.shapeData.iteritems():
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
				else: 
					print("-")
			exp_offset += sd.reactivity.shape[0]
			n = 0

		collection = PatchCollection(patches,match_original=True)
		plt.add_collection(collection)
		plt.axis([i,0,0,exp_offset])


	def reactivity_plot(self,plt):
		x = range(len(self.shapeData[self.select_reactivity].sequence))
		color_bars = empty(len(x),dtype="S10")
		react_color = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
		text_cooloff = 0
		for n,base in enumerate(self.shapeData[self.select_reactivity].sequence):
			flex_class = int(self.shapeData[self.select_reactivity].mean_reactivity_class[n])
			color_bars[n] = str(react_color[flex_class])
			if ((flex_class >= 3 or n%20==0) and text_cooloff<=0 ):
				text_cooloff = 4
				plt.text(n-0.5,self.shapeData[self.select_reactivity].mean_flex[n],str(base)+str(n+1),color=str(react_color[flex_class]))
			text_cooloff -= 1

		plt.bar(x,self.shapeData[self.select_reactivity].mean_flex,color=color_bars,ec="none")
		plt.set_xlim([n,0])
		for c,pos_y in enumerate((0.2,0.4,0.6,0.8,1)):
			plt.axhline(y=pos_y,color=react_color[c])
			plt.text(n-1,pos_y,str(pos_y),color=react_color[c])
		

	def sdreactivity_plot(self,plt):
		x = range(len(self.shapeData[self.select_sdreactivity].sequence))
		color_bars = empty(len(x),dtype="S10")
		react_color = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#AA0078"] 
		text_cooloff = 0
		for n,base in enumerate(self.shapeData[self.select_sdreactivity].sequence):
			color = int(self.shapeData[self.select_sdreactivity].variability_class[n])
			if color > 5:
				color = 5
			color_bars[n] = str(react_color[color])
			if ((color >= 3 or n%20==0) and text_cooloff<=0 ):
				text_cooloff = 4
				plt.text(n-0.5,self.shapeData[self.select_sdreactivity].var[n],str(base)+str(n+1),color=str(react_color[color]))
			text_cooloff -= 1

		plt.bar(x,self.shapeData[self.select_sdreactivity].var,color=color_bars,ec="none")
		plt.set_xlim([n,0])
		for i in range(5):
			pos_y=self.shapeData[self.select_sdreactivity].varmean+i*self.shapeData[self.select_sdreactivity].sd_of_var
			plt.axhline(y=pos_y,color=react_color[i])
			plt.text(n-1,pos_y,"Mean + "+str(i)+" times the standard deviation",color=react_color[i])
		
 


class ModelHandler(Handler):
    panel = Instance(ModelPanel)
    
    def __init__(self,panel):
        self.panel = panel
