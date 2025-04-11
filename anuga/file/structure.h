
#include <stdint.h>

/*** new RWG DEFINITION FOR TIDE GAUGE OUTPUT ***/
struct tgsrwg
   {
   float geolat;
   float geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
   float mcolat;
   float mcolon;
   int32_t ig;
   int32_t ilon;     /* grid point location */
   int32_t ilat;
   float z;             /* water depth at this location */
   float center_lat, center_lon;	/**** center of this array *****/
   float offset,az,baz;			/* for arrays this is the distance in km from center of array */
   float dt;		/* sampling rate for this site */
   int32_t nt;		/* number of points */
   char id[16];         /* identifier */
   };

/*** DEFINITION FOR TIDE GAUGE OUTPUT ***/
struct tgs
   {
   int32_t ista;
   float geolat,geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
   int32_t   ilat,ilon;     /* grid point location */

   float z;             /* water depth at this location */

   float center_lat, center_lon;	/**** center of this array *****/
   float offset,az,baz;			/* for arrays this is the distance in km from center of array */

   float dt;		/* sampling rate for this site */
   int32_t nt;		/* number of points */

   char id[16];         /* identifier */
   };


