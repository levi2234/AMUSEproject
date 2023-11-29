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
    
        
    def plot(self, save=True,c="blue", close=True,**kwargs):
        for systemobject in self.system: 
            xaxisval = self.dynamic_accessor(systemobject, self.xaxis).value_in(units.km)
            yaxisval = self.dynamic_accessor(systemobject, self.yaxis).value_in(units.km)
            
            try: 
                objectlabel = self.dynamic_accessor(systemobject, "name")
            except: 
                objectlabel = None
            plt.scatter(xaxisval,yaxisval,label = objectlabel,c=c,**kwargs)
            
        if close==False:
            return  
            
        if self.legend:
            plt.legend()
        if self.grid:
            plt.grid()
        
 
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if isinstance(save, str):
            plt.savefig(save)
        

            
        
        
        plt.show()
        
        
    def __init__(self,system, xaxis = "x", yaxis = "y", xlabel=None, ylabel=None, legend=False, grid=False): #init is placed last so we can call the plot function in init
        self.system = system
        self.xaxis = xaxis
        self.yaxis = yaxis

        
        
        if xlabel == None:
            self.xlabel = xaxis
        else: 
            self.xlabel = xlabel
            
        if ylabel == None:
            self.ylabel = yaxis
            
        else: 
            self.ylabel = ylabel

        self.legend = legend
        self.grid = grid
 
        self.plot()

