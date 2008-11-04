# setup.py
#run with "python exesetup.py py2exe"
from distutils.core import setup
import py2exe
      
setup(console=["pmesh.py"],
      data_files=[("icons",
                   ["icons/addVertex.gif",
                    "icons/autoSeg.gif",
                    "icons/edit.gif",
                    "icons/hole.gif",
                    "icons/joinVer.gif",
                    "icons/meshGen.gif",
                    "icons/pointer.gif",
                    "icons/default.gif",
                    "icons/region.gif",
                    "icons/segment.gif",
                    "icons/sep.gif",
                    "icons/vertex.gif",
                    "icons/zoom0.5.gif",
                    "icons/zoom1.0.gif",
                    "icons/zoom2.gif",
                    "icons/zoomToMesh.gif"]),
                  (".",["./hull.exe"])])
