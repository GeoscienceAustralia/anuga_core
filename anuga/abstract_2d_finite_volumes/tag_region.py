"""region.py - Classes for implementing region conditions

NOTE: I've only set it up for testing X values of constants and functions.
Although it should work for vectors/arrays of values, that implies
knowing info about the actual triangles, and that is not how the user should
operate.
"""
# FIXME (DSG-DSG) add better comments

from builtins import object
import numpy as num


class Tag_region(object):
    """Base class for modifying quantities based on a region.
    """   

    def __init__(self, location='vertices'):
        self.location = location
        
    def __call__(self, tag, elements, domain):
        msg = 'Generic class Boundary must be subclassed'
        raise Exception(msg)


    def build_indices(self, elements, domain):
        """
        Return a list of triangle_id or vertex_id, depending on the location
        """
        if self.location == 'unique vertices':
            return domain.get_unique_vertices(elements)
        else:
            return elements
                
class Set_tag_region(Tag_region):
    
    def __init__(self, tag, quantity, X, location='vertices'):
        """ 
        name: Name of quantity
        X: const or function 
        location: Where values are to be stored.
        Permissible options are: vertices, centroid
        """
        
        Tag_region.__init__(self)
        self.tag = tag 
        self.quantity = quantity
        self.location = location
        self.X = X

    def __repr__(self):
        pass
    
    def __call__(self, tag, elements, domain):
        """
        """    
        if tag == self.tag:
            domain.set_quantity(self.quantity,
                                self.X,
                                location=self.location,
                                indices=self.build_indices(elements, domain))

        
class Add_value_to_region(Tag_region):
    """
    Will add a value to the current quantity value.
    
    quantity = initial_quantity + X.
    
    This method does not work with functions
    
    Use average to add X to a mean average of the quantity values.
    eg if applying this to elevation this will give a flat roof.
    """
    
    def __init__(self, tag, quantity, X, location='vertices', initial_quantity=None, average=False):
        #I have to get this going!
        #Region.__init__(self)
        self.tag = tag 
        self.quantity_answer = quantity
        self.location = location
        self.X = X
        self.average = average
        if initial_quantity is None:
            self.quantity_initial_value = quantity
        else:
            self.quantity_initial_value = initial_quantity
        if callable(X):
            raise Exception('This class does not work with functions')

    def __repr__(self):
        pass
    
    def __call__(self, tag, elements, domain):
        """
        """    
        if tag == self.tag:
            #new_values = domain.get_quantity(self.quantity_initial_value,
            #              indices=self.build_indices(elements, domain),
            #              location=self.location) + self.X
            Q = domain.get_quantity(self.quantity_initial_value)
            if self.average is True:
                # Average the points, and then the verts in the averaged point.  
                values = Q.get_values(indices=self.build_indices(elements, domain),
                                      location=self.location)
                av = num.average(values)
                if self.location == "vertices":
                    av = num.average(av)
                new_values = av + self.X    
            else:
                new_values = Q.get_values(indices=self.build_indices(elements, domain),
                                      location=self.location) + self.X    
            domain.set_quantity(self.quantity_answer, new_values,
                                indices=self.build_indices(elements, domain),
                                location=self.location)

class Add_quantities(Tag_region):
    """
    Will add a quantity to the current quantity value.
    """
    
    def __init__(self, tag, quantity_answer, adding_quantity, location='vertices'):
        #I have to get this going!
        #Region.__init__(self)
        self.tag = tag 
        self.quantity_answer = quantity_answer
        self.adding_quantity = adding_quantity
        self.location = location

    def __repr__(self):
        pass
    
    def __call__(self, tag, elements, domain):
        """
        """    
        if tag == self.tag:

            #new_values = domain.get_quantity(self.quantity_answer,
            #              indices=self.build_indices(elements, domain),
            #              location=self.location) \
            #              + domain.get_quantity(self.adding_quantity,
            #              indices=self.build_indices(elements, domain),
            #              location=self.location)



            indices = self.build_indices(elements, domain)
            location = self.location
            Q1 = domain.get_quantity(self.quantity_answer)
            Q2 = domain.get_quantity(self.adding_quantity)            
            new_values = Q1.get_values(indices=indices, location=location) +\
                         Q2.get_values(indices=indices, location=location)

            
            domain.set_quantity(self.quantity_answer, new_values,
                                indices=self.build_indices(elements, domain),
                                location=self.location)


class Stage_no_less_than_elevation(Tag_region):
    """
    Will set the stage to not be less than the elevation.
    This would be good, but it's not region dependent.
    Wait for it to become a default for pyvolution.
    """
    
    def __init__(self):
        pass
