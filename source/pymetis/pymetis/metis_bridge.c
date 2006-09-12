/* Python.h and metis.h don't get along, so this is an ugly
 * hack to make it work.
 */

//#include <metis.h>

#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <rename.h>
#include <proto.h>


void bridge_partMeshNodal(int * ne, int * nn, idxtype * elmnts, int * etype, int * numflag, int * nparts, int * edgecut, idxtype * epart, idxtype * npart){
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
