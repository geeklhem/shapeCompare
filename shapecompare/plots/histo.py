from numpy import empty
palette = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
bbox_props = dict(boxstyle="square", fc="w", ec="0.5", alpha=0.9)


def reactivity_mean(plt,data):
        x = range(len(data.sequence))
        color_bars = empty(len(x),dtype="S10")
        text_cooloff = 0
        for n,base in enumerate(data.sequence):
                flex_class = int(data.mean_reactivity_class[n])
                color_bars[n] = str(palette[flex_class])
                if ((flex_class >= 3 or n%20==0) and text_cooloff<=0 ):
                        text_cooloff = 2
                        if data.mean_flex[n]>0:
                            y = data.mean_flex[n]+0.1
                        else:
                            y = data.mean_flex[n]-0.1
                        plt.text(n+1,
                                 y,
                                 str(base)+str(n+1),
                                 color=str(palette[flex_class]),
                                 rotation="vertical")
                text_cooloff -= 1

        plt.bar(x,data.mean_flex,yerr=data.var,ecolor="grey",color=color_bars,ec="none")
        plt.set_xlim([n,0])
        plt.text(0.03,0.95,data.name+" mean reactivity",
                 transform=plt.transAxes,
                 bbox=bbox_props)
        

        for c,pos_y in enumerate((0.2,0.4,0.6,0.8,1)):
                plt.axhline(y=pos_y,color=palette[c])
                plt.text(n-1,pos_y,str(pos_y),color=palette[c])
        

def reactivity_uniq(plt,data,exp):
        x = range(len(data.sequence))
        color_bars = empty(len(x),dtype="S10")
        text_cooloff = 0
        for n,base in enumerate(data.sequence):
                flex_class = int(data.reactivities_class[exp,n])
                color_bars[n] = str(palette[flex_class])
                if ((flex_class >= 3 or n%20==0) and text_cooloff<=0 ):
                        text_cooloff = 2
                        if data.mean_flex[n]>0:
                            y = data.reactivity[exp,n]+0.1
                        else:
                            y = data.reactivity[exp,n]-0.1
                        plt.text(n+1,
                                 y,
                                 str(base)+str(n+1),
                                 color=str(palette[flex_class]),
                                 rotation="vertical")
                text_cooloff -= 1

        plt.bar(x,data.reactivity[exp],yerr=data.var,ecolor="grey",color=color_bars,ec="none")
        plt.set_xlim([n,0])
        plt.text(0.03,0.95,"["+data.name+"]:"+data.names[exp]+" mean reactivity",
                 transform=plt.transAxes,
                 bbox=bbox_props)
        

        for c,pos_y in enumerate((0.2,0.4,0.6,0.8,1)):
                plt.axhline(y=pos_y,color=palette[c])
                plt.text(n-1,pos_y,str(pos_y),color=palette[c])


def scatter(plt,data):
        x = range(len(data.sequence))
        colors = []       
        for flex_class in data.mean_reactivity_class: 
            colors.append(str(palette[int(flex_class)]))
                
        z = zip(data.mean_flex,data.var,colors)
        z.sort()
        mean_flex, var, colors = zip(*z)
        plt.errorbar(x,mean_flex,yerr=var,ecolor="grey",ls='none')
        plt.scatter(x,mean_flex,c=colors,marker="o")
        plt.set_xlim([x[-1],0])
        plt.set_ylim(0,1)
        plt.text(0.03,0.95,data.name+" mean reactivity categories",
                 transform=plt.transAxes,
                 bbox=bbox_props)

        for c,pos_y in enumerate((0.2,0.4,0.6,0.8,1)):
                plt.axhline(y=pos_y,color=palette[c])
                plt.text(1,pos_y,str(pos_y),color=palette[c])
        

