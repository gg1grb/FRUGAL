c--------------------------"transform.h"------------------------------
c  contains the temporary arrays used in transforming data to 
c  model grid
c--------------------------------------------------------------------
      real datain(-180:180,-90:90),datainv(-180:180,-90:90,2)
#if defined daily || defined nordic_tracer || defined nordic_heateof
      real dataday(192,-46:47),datadayv(192,-46:47,2)
#endif
     &    ,datamodel(imt,jmt),datamodelv(imt,jmt,2)
      common datain,datainv,datamodel,datamodelv
#if defined daily || defined nordic_tracer || defined nordic_heateof
     &    ,dataday,datadayv
#endif
