
#-------------------------------------------------------------
if __name__ == "__main__":

        import tempfile
        import os
    
        from anuga.shallow_water import Domain, Reflective_boundary, \
                            Dirichlet_boundary
       
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        
        # Create a scenario outline.
        polygon = [[0,0],[100,0],[100,100],[0,100]]
        
        boundary_tags = {'wall':[0,1,3],'wave':[2]}
        
        inner1_polygon = [[10,10],[20,10],[20,20],[10,20]]
        

        inner2_polygon = [[30,30],[40,30],[40,40],[30,40]]
        
        
        max_area = 1
        interior_regions = [(inner1_polygon, 5),(inner2_polygon, 10)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions)

        # Create an ungenerate file            
        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write("         1       ??      ??\n\
       90.0       90.0\n\
       81.0       90.0\n\
       81.0       81.0\n\
       90.0       81.0\n\
       90.0       90.0\n\
END\n\
         2      ?? ??\n\
       10.0       80.0\n\
       10.0       90.0\n\
       20.0       90.0\n\
       10.0       80.0\n\
END\n\
END\n")
        file.close() 

        # import the ungenerate file
        m.import_ungenerate_file(fileName) 
        os.remove(fileName)
	
        m.generate_mesh(maximum_triangle_area=max_area,verbose=False)
        mesh_filename = "mesh.tsh"
        m.export_mesh_file(mesh_filename)

	# Run a simulation on the mesh
        domain = Domain(mesh_filename, use_cache = False)
        
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([3,0,0]) 
        domain.set_boundary( {'wall': Br, 'wave': Bd} )
        yieldstep = 0.1
        finaltime = 20
        for t in domain.evolve(yieldstep, finaltime):    
            domain.write_time()

            
