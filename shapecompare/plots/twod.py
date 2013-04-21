palette = ['#66CCFF','#33FF33','#FFFF33','#FF6600',"#FF0000","#CCCCCC","#AA0078"] 
bbox_props = dict(boxstyle="square", fc="w", ec="0.5", alpha=0.9)

def mean(plt,model,data):
    app = 0
    for n,(xpos,ypos) in enumerate(model.position):
        color = int(data.mean_reactivity_class[n])
        if model.app[n]:
            app = "bold"
        else:
            app = "normal"
        plt.text(xpos,ypos,
                 str(data.sequence[n]),
                 backgroundcolor=palette[color],
                 weight=app)

    plt.axis([0,1,1,0])
    plt.axis('equal')
    plt.text(0.03,0.95,data.name+" mean reactivity on "+model.name+" model.",
             transform=plt.transAxes,
             bbox=bbox_props)


def uniq(plt,model,data,exp):
    app = 0
    for n,(xpos,ypos) in enumerate(model.position):
        color = int(data.reactivities_class[exp,n])
        if model.app[n]:
            app = "bold"
        else:
            app = "normal"
        plt.text(xpos,ypos,
                 str(data.sequence[n]),
                 backgroundcolor=palette[color],
                 weight=app)

    plt.axis([0,1,1,0])
    plt.axis('equal')
    plt.text(0.03,0.95,data.names[exp]+"reactivity on "+model.name+" model.",
             transform=plt.transAxes,
             bbox=bbox_props)
