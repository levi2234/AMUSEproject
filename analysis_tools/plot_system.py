import matplotlib.pyplot as plt
import amuse.units as u
import amuse.units as units
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

    def __init__(self,system, xaxis = "x", yaxis = "y", xlabel=None, ylabel=None, legend=False, grid=False):
        self.system = system
        self.xaxis = xaxis
        self.yaxis = yaxis
        
        if xlabel == None:
            xlabel = xaxis
        else: 
            xlabel = xlabel
            
        if ylabel == None:
            ylabel = yaxis
            
        else: 
            ylabel = ylabel

        self.legend = legend

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
    
        
    def plot(self,**kwargs):
        for systemobject in self.system: 
            xaxisval = self.dynamic_accessor(systemobject, self.xaxis).value_in(units.AU)
            yaxisval = self.dynamic_accessor(systemobject, self.yaxis).value_in(units.AU)
            try: 
                objectlabel = self.dynamic_accessor(systemobject, "name")
            except: 
                objectlabel = None
            plt.scatter(xaxisval,yaxisval,label = objectlabel,**kwargs)
            
        if self.legend:
            plt.legend()
        if self.grid:
            plt.grid()
            
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.show()