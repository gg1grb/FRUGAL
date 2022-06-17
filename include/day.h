c---------------------------------------
c  this code defines the 2922 days of the daily forcing fields
c---------------------------------------
      CHARACTER*1 SEG(4),UNIT
      INTEGER COUNT,Ldd,Idd,Jdd,Kdd
      COUNT=0
      DO Ldd=0,2
      IF(Ldd.eq.0)THEN
      SEG(1)='0'
      ELSEif(Ldd.eq.1)then
      SEG(1)='1'
      else
      seg(1)='2'
      ENDIF
        DO Idd=0,9
        SEG(2)=UNIT(Idd)
          DO Jdd=0,9
          SEG(3)=UNIT(Jdd)
            DO Kdd=0,9
            COUNT=COUNT + 1
            IF(COUNT.eq.2923) goto 388
            SEG(4)=UNIT(Kdd)
            iDAY(COUNT)=SEG(1)//SEG(2)//SEG(3)//SEG(4)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
388   continue
