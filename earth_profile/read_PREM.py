from amuse.units import units,constants
import numpy as np
import csv

class CSVReader:
    """
    A class for reading the Preliminary Earth
    Reference Module
    Attributes:
    ----------
    filename: it should be PREM.csv
    """
    
    def __init__(self, filename):
        self.filename = filename
        self.data_dict = {'radius_profile': [],'density_profile': [],
                          'radius':6371000|units.m, 'mass': 5.972e27 | units.g}

    def read_and_convert(self):
        """
        Reads the PREM csv file and stores the density and radial profiles
        in a dictionary with AMUSE units
        """
        with open(self.filename, 'r') as file:
            reader = csv.reader(file)
            # Skip the first two rows
            next(reader)
            next(reader)

            # Extract and convert the first column to AMUSE units of meters
            for row in reader:
                if row:
                    radius_in_km = float(row[0])
                    density_in_gcm3 = float(row[2])
                    self.data_dict['radius_profile'].append(1000*radius_in_km)
                    self.data_dict['density_profile'].append(density_in_gcm3)
            self.data_dict['radius_profile']=self.data_dict['radius_profile']|units.m
            self.data_dict['density_profile']=self.data_dict['density_profile']|units.g * units.cm**(-3)
    def get_data(self):
        """
        Returns a dictionary with the Earth's radial and density profiles

        Returns:
            Dictionary with the following keys:
            - 'radius_profile': radial profile in AMUSE units of meters
            - 'density_profile': density profile in AMUSE units of g*cm**-3
            - 'mass': total mass of the Earth in AMUSE units of g
            - 'radius': total radius of the Earth in AMUSE units of m
        """
        return self.data_dict
 
