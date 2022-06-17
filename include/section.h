c------------------------"section.h"----------------------------
c contains variables for lat-long section location
c--------------------------------------------------------------

      integer sec(imt*2+jmt*2,2)
      integer ninsec
      common /sectiondata/ sec,ninsec
