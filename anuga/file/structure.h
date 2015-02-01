
/*** new RWG DEFINITION FOR TIDE GAUGE OUTPUT ***/
struct tgsrwg
   {
   float geolat;
   float geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
   float mcolat;
   float mcolon;
   int ig;
   int ilon;     /* grid point location */
   int ilat;
   float z;             /* water depth at this location */
   float center_lat, center_lon;	/**** center of this array *****/
   float offset,az,baz;			/* for arrays this is the distance in km from center of array */
   float dt;		/* sampling rate for this site */
   int nt;		/* number of points */
   char id[16];         /* identifier */
   };

/*** DEFINITION FOR TIDE GAUGE OUTPUT ***/
struct tgs
   {
   int ista;
   float geolat,geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
   int   ilat,ilon;     /* grid point location */

   float z;             /* water depth at this location */

   float center_lat, center_lon;	/**** center of this array *****/
   float offset,az,baz;			/* for arrays this is the distance in km from center of array */

   float dt;		/* sampling rate for this site */
   int nt;		/* number of points */

   char id[16];         /* identifier */
   };


