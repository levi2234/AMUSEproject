import matplotlib.pyplot as plt
from amuse.units import units  
#function for plotting an amuse particle system

class systemplotter:
    """
    A class for creating scatter plots of system objects.

    Attributes:
    -----------
    system : Amuse ParticleSet
        A list of system objects to be plotted (e.g. planets, stars, moons, etc.) 
    xaxis : str
        The attribute name to be used for the x-axis.
    yaxis : str
        The attribute name to be used for the y-axis.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.
    legend : bool
        Whether or not to display a legend.
    grid : bool
        Whether or not to display a grid.
    """


    def dynamic_accessor(self, obj, attribute_name):
        """
        Returns the value of the specified attribute of the given object, if it exists.

        Args:
            obj: The object to retrieve the attribute from.
            attribute_name: The name of the attribute to retrieve.

        Returns:
            The value of the specified attribute.

        Raises:
            AttributeError: If the specified attribute does not exist on the object.
        """
        if hasattr(obj, attribute_name):
            return getattr(obj, attribute_name)
        else:
            raise AttributeError(f"'{type(obj).__name__}' object has no attribute '{attribute_name}'")
    
        

    def plot(self, savefig=False, **kwargs):
        fig, ax = plt.subplots()

        for systemobject in self.system: 
            xaxisval = self.dynamic_accessor(systemobject, self.xaxis).value_in(units.m)
            yaxisval = self.dynamic_accessor(systemobject, self.yaxis).value_in(units.m)
            try: 
                objectlabel = self.dynamic_accessor(systemobject, "name")
            except: 
                objectlabel = None
                
            if objectlabel == "planet":
                radius = self.dynamic_accessor(systemobject, "radius").value_in(units.m)
                circle = plt.Circle((0, 0), radius, color='green',alpha= 0.2, fill=True)
                ax.add_artist(circle)
            else: 
                ax.scatter(xaxisval, yaxisval, label=objectlabel, **kwargs)
            
        if self.legend:
            ax.legend()
        if self.grid:
            ax.grid()
            
        #ax.set_aspect('equal')
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        if savefig:
            plt.savefig(savefig)
        else:
            plt.show()
        
    def __init__(self, system, xaxis="x", yaxis="y", xlabel=None, ylabel=None, legend=False, grid=False, **kwargs):
        self.system = system
        self.xaxis = xaxis
        self.yaxis = yaxis
        
        
        if xlabel is None:
            self.xlabel = xaxis
        else: 
            self.xlabel = xlabel
            
        if ylabel is None:
            self.ylabel = yaxis
        else: 
            self.ylabel = ylabel

        self.legend = legend
        self.grid = grid
        
        self.plot()
