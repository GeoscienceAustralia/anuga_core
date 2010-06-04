 
##
# @brief Convert points to a polygon (?)
# @param points_file The points file.
# @param minimum_triangle_angle ??
# @return 
def load_pts_as_polygon(points_file, minimum_triangle_angle=3.0):
    """
    WARNING: This function is not fully working.

    Function to return a polygon returned from alpha shape, given a points file.

    WARNING: Alpha shape returns multiple polygons, but this function only
             returns one polygon.
    """

    from anuga.pmesh.mesh import Mesh, importMeshFromFile
    from anuga.shallow_water.shallow_water_domain import Domain
    from anuga.pmesh.mesh_interface import create_mesh_from_regions

    mesh = importMeshFromFile(points_file)
    mesh.auto_segment()
    mesh.exportASCIIsegmentoutlinefile("outline.tsh")
    mesh2 = importMeshFromFile("outline.tsh")
    mesh2.generate_mesh(maximum_triangle_area=1000000000,
                        minimum_triangle_angle=minimum_triangle_angle,
                        verbose=False)
    mesh2.export_mesh_file('outline_meshed.tsh')
    domain = Domain("outline_meshed.tsh", use_cache = False)
    polygon =  domain.get_boundary_polygon()
    return polygon
