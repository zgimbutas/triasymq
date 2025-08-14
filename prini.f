cc Copyright (C) 2014: Vladimir Rokhlin
cc
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory).
cc
cc SPDX-License-Identifier: BSD-3-Clause-Modification
cc
cc Modified August 13, 2025: Zydrunas Gimbutas
cc - remove unnecessary save statements
cc - by default, initialize printing units to 0
cc - remove fileflush, mach_zero and ztime
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Printing routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
        SUBROUTINE PRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
        INTEGER IA(1)
        INTEGER *2 IA2(1)
        LOGICAL *1 LA(1)
        data IP/0/,IQ/0/
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
C
C
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C
C
C
C
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C
C
C
C
        ENTRY PRIN2_long(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1450)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1450)(A2(J),J=1,N)
 1450 FORMAT(2(2X,E22.16))
        RETURN
C
C
C
C
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C
C
C
C
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C
C
C
C
        ENTRY PRINF1(MES,IA1,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA1(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA1(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINL1(MES,LA,N)
        CALL MESSPR(MES,ip,iq)
 2200 format(20L3)
c
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2200)(lA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2200)(lA(J),J=1,N)
c
        RETURN
        END
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I1=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
c
c
c
c
c
        subroutine msgmerge(a,b,c)
        character *1 a(*),b(*),c(*),ast
        data ast/'*'/
c
        do 1200 i=1,1000
c
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)
        iadd=i
 1200 continue
c
 1400 continue
c
        do 1800 i=1,1000
c
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c
