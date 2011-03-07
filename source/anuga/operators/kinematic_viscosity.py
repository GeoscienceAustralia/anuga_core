from anuga import Domain
#from anuga.utilities.sparse import Sparse, Sparse_CSR
from sparse import Sparse, Sparse_CSR #the new import
from anuga.utilities.cg_solve import conjugate_gradient
import anuga.abstract_2d_finite_volumes.neighbour_mesh as neighbour_mesh
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Dirichlet_boundary
import numpy as num
import kinematic_viscosity_ext
import anuga.utilities.log as log

class Kinematic_Viscosity_Operator:

	def __init__(self, domain, triangle_areas=True, verbose=False):
		if verbose: log.critical('Kinematic Viscosity: Beginning Initialisation')
		#Expose the domain attributes
		self.domain = domain
		self.mesh = domain.mesh
		self.boundary = domain.boundary
		self.n = len(self.domain)
		self.dt = 1e-6 #Need to set to domain.timestep
		self.boundary_len = len(domain.boundary)
		self.tot_len = self.n+self.boundary_len
		self.boundary_enum = self.enumerate_boundary()
		self.verbose = verbose
		
		#Geometric Information
		if verbose: log.critical('Kinematic Viscosity: Building geometric structure')
		self.geo_structure_indices = num.zeros((self.n,3),num.int)
		self.geo_structure_values = num.zeros((self.n,3),num.float)
		kinematic_viscosity_ext.build_geo_structure(self,self.n,self.tot_len)
		self.apply_triangle_areas = triangle_areas
		self.triangle_areas = Sparse(self.n,self.n)
		for i in range(self.n):
			self.triangle_areas[i,i] = 1.0 / self.mesh.areas[i]
		self.triangle_areas = Sparse_CSR(self.triangle_areas)
		
		#Application Information
		self.qty_considered = 1 #1 or 2 (uh or vh respectively)
		self.operator_matrix = Sparse(self.n,self.tot_len)
		self.operator_data = num.zeros((4*self.n,),num.float) #Sparse_CSR.data
		self.operator_colind = num.zeros((4*self.n,),num.int) #Sparse_CSR.colind
		self.operator_rowptr = 4*num.indices((self.n+1,))[0,:] #Sparse_CSR.rowptr (4 entries in every row, we know this already) = [0,4,8,...,4*n]
		self.stage_heights = num.zeros((self.n,1),num.float)
		self.stage_heights_scaling = Sparse(self.n,self.n)
		self.boundary_vector = num.zeros((self.n,),num.float)
		
		self.parabolic_solve = False #Are we doing a parabolic solve at the moment?
		
		if verbose: log.critical('Kinematic Viscosity: Initialisation Done')
    
	# Enumerate the boundary conditions in some way
	def enumerate_boundary(self):
		#Just enumerate by the triangle number, then edge number
		enumeration = {}
		enum = 0
		for i in range(self.n):
			for edge in range(3):
				if self.mesh.neighbours[i,edge] == -1:
					enumeration[(i,edge)] = enum
					enum += 1
		return enumeration
	
	def set_qty_considered(self,qty):
		if qty == 1 or qty == 'u':
			self.qty_considered = 1
		elif qty == 2 or qty == 'v':
			self.qty_considered = 2
		else: #Raise an exception
			msg = "Incorrect input qty"
			assert 0==1, msg
	
	def apply_stage_heights(self,h):
		msg = "(KV_Operator.apply_stage_heights) h vector has incorrect length"
		assert h.size == self.n, msg
		
		self.operator_matrix = Sparse(self.n,self.tot_len)
		
		#Evaluate the boundary stage heights - use generic_domain.update_boundary? (check enumeration)
		boundary_heights = num.zeros((self.boundary_len,1),num.float)
		for key in self.boundary_enum.keys():
			boundary_heights[self.boundary_enum[key]] = self.boundary[key].evaluate()[0]
		
		kinematic_viscosity_ext.build_operator_matrix(self,self.n,self.tot_len,h,boundary_heights)
		self.operator_matrix = Sparse_CSR(None,self.operator_data,self.operator_colind,self.operator_rowptr,self.n,self.tot_len)
		
		self.stage_heights = h
		#Set up the scaling matrix
		data = h
		num.putmask(data, data!=0, 1/data) #take the reciprocal of each entry unless it is zero
		self.stage_heights_scaling = Sparse_CSR(None,data,num.arange(self.n),num.append(num.arange(self.n),self.n),self.n,self.n)
		
		self.build_boundary_vector()
		
		return self.operator_matrix
	
	def build_boundary_vector(self):
        #Require stage heights applied
		X = num.zeros((self.tot_len,2),num.float)
		for key in self.boundary_enum.keys():
			quantities = self.boundary[key].evaluate()
			h_i = quantities[0]
			if h_i == 0.0:
				X[self.n + self.boundary_enum[key],:] = num.array([0.0, 0.0])
			else:
				X[self.n + self.boundary_enum[key],:] = quantities[1:] / h_i
		self.boundary_vector = self.operator_matrix * X
		#Tidy up
		if self.apply_triangle_areas:
			self.boundary_vector = self.triangle_areas * self.boundary_vector

		return self.boundary_vector
	
	def elliptic_multiply(self,V,qty_considered=None,include_boundary=True):
		msg = "(KV_Operator.apply_vector) V vector has incorrect dimensions"
		assert V.shape == (self.n,) or V.shape == (self.n,1), msg
		
		if qty_considered!=None:
			print "Quantity Considered changed!"
			self.set_qty_considered(qty_considered)
		
		D = num.zeros((self.n,),num.float)
		#Change (uh) to (u), (vh) to (v)
		X = self.stage_heights_scaling * V
		
		if self.apply_triangle_areas:
			X = self.triangle_areas * X
		
		X = num.append(X,num.zeros((self.boundary_len,)),axis=0)
		D = self.operator_matrix * X
		
		if include_boundary:
			D += self.boundary_vector[:,self.qty_considered-1]
		
		D = self.stage_heights * D #make sure we return (uh) not (u), for the solver's benefit
		
		return num.array(D).reshape(self.n,)
	
	#parabolic_multiply(V) = identity - dt*elliptic_multiply(V)
	#We do not need include boundary, as we will only use this method for solving (include_boundary = False)
	#Here V is already either (u) or (v), not (uh) or (vh)
	def parabolic_multiply(self,V):
		msg = "(KV_Operator.parabolic_multiply) V vector has incorrect dimensions"
		assert V.shape == (self.n,) or V.shape == (self.n,1), msg
		
		D = V #The return value
		X = V #a temporary array
		
		#Apply triangle areas
		if self.apply_triangle_areas:
			X = self.triangle_areas * X
		
		#Multiply out
		X = num.append(X,num.zeros((self.boundary_len,)),axis=0)
		D = D - self.dt * (self.operator_matrix * X)
		
		return num.array(D).reshape(self.n,)
		
	def __mul__(self,other):
		try:
			B = num.array(other)
		except:
			msg = "Trying to multiply the Kinematic Viscosity Operator against a non-numeric type"
			raise msg
		
		if len(B.shape) == 0:
			#Scalar
			R = B*self
		elif len(B.shape) == 1 or (len(B.shape) == 2 and B.shape[1]==1):
			#Vector
			if self.parabolic_solve:
				R = self.parabolic_multiply(other)
			else:
				#include_boundary=False is this is *only* used for cg_solve()
				R = self.elliptic_multiply(other,include_boundary=False)
		else:
			raise ValueError, 'Dimension too high: d=%d' %len(B.shape)
		return R
    
	def __rmul__(self, other):
		#Right multiply with scalar
		try:
			other = float(other)
		except:
			msg = 'Sparse matrix can only "right-multiply" onto a scalar'
			raise TypeError, msg
		else:
			new = self.operator_matrix * new
		return new
	
	def cg_solve(self,B,qty_to_consider=None):
		if len(B.shape) == 1:
			return self.cg_solve_vector(B,qty_to_consider)
		elif len(B.shape) == 2:
			return self.cg_solve_matrix(B)
		else:
			raise ValueError, 'Dimension too high: d=%d' %len(B.shape)
	
	def cg_solve_matrix(self,B):
		assert B.shape[1] < 3, "Input matrix has too many columns (max 2)"
		
		X = num.zeros(B.shape,num.float)
		for i in range(B.shape[1]):
			X[:,i] = self.cg_solve_vector(B[:,i],i+1) #assuming B columns are (uh) and (vh)
		return X
    
	def cg_solve_vector(self,b,qty_to_consider=None):
		if not qty_to_consider==None:
			self.set_qty_considered(qty_to_consider)
		x = conjugate_gradient(self,b - self.boundary_vector[:,self.qty_considered-1]) #Call the ANUGA conjugate gradient utility
		return x
	
	#Solve the parabolic equation to find u^(n+1), where u = u^n
	#Here B = [uh vh] where uh,vh = n*1 column vectors
	def parabolic_solver(self,B):
		assert B.shape == (self.n,2), "Input matrix has incorrect dimensions"
		
		next_B = num.zeros((self.n,2),num.float)
		#The equation is self.parabolic_multiply*u^(n+1) = u^n + (dt)*self.boundary_vector
		self.parabolic_solve = True
		#Change (uh) and (vh) to (u) and (v)
		b = self.stage_heights_scaling * B

		next_B = conjugate_gradient(self, b + (self.dt * self.boundary_vector), iprint=1)
		
		#Return (u) and (v) back to (uh) and (vh)
		next_B = self.stage_heights_scaling * next_B
		
		self.parabolic_solve = False
		
		return next_B
