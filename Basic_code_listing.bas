'In Procsides it does not matter which way round the vertex coords 
'are used,i.e., (x1-x3) or (x3-x1) because the square is always positive
OPTION EXPLICIT
DECLARE SUB ProcForward ()
DECLARE SUB ProcTryLoops ()
DECLARE SUB ProcVertexCoords ()
DECLARE SUB ProcCircumscribed ()
DIM SHARED AS UINTEGER knt
DIM SHARED AS DOUBLE s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 'sides to calculate
DIM SHARED AS DOUBLE x,y,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9 'vertices to calculate
DIM SHARED AS DOUBLE a,b,c,d,e,f,h,aC,aB,tbl,tblx,arm,armx'vari constants
DIM SHARED AS DOUBLE astart,bstart,cstart,dstart,estart,fstart,hstart,armxstart,tblxstart,range,ink
DIM SHARED AS DOUBLE armxend,tblxend,aend,bend,cend,dend,eend,fend,hend
DIM SHARED AS DOUBLE tinc0,tinc1,tinc2,tinc3,tinc4,tinc5,tinc6,tinc7,tinc8,tinc9,nearest'constants
DIM SHARED AS DOUBLE ainc0,ainc1,ainc2,ainc3,ainc4,ainc5,ainc6,ainc7,ainc8,ainc9',lo,hi'constants
tinc1= 199.764258931549    :ainc1= 150.699018242643   '<<<<<<<<<<<<table and arm screw incs
tinc2= 137.477763809118    :ainc2= 181.699831314401
tinc3=  54.74804692617     :ainc3= 172.900372581192
tinc4=   1.112002513305    :ainc4= 128.951849742162
tblx=  249.745657625133    :ainc5=  72.29493974962
tinc6=  39.569822405649    :ainc6=  24.40537191706
tinc7=  98.568865748128    :armx=  222.511658538965
tinc8= 157.817063265319    :ainc8=   6.477230466509
tinc9= 202.011409140245    :ainc9=  42.444564873512
tinc0= 219.562952509625    :ainc0=  96.684983717207

a=745:b=745:c=745:d=694:e=692:f=583.125:h=120
tblx=249.745657625133:armx=222.511658538965
range=1.5
ink=.5
'lo=165.995:hi=166.005
nearest=1 
astart=a   -range:   aend=a   +range
bstart=b   -range:   bend=b   +range
cstart=c   -range:   cend=c   +range
dstart=d   -range:   dend=d   +range
estart=e   -range:   eend=e   +range
fstart=f   -range:   fend=f   +range
hstart=h   -range:   hend=h   +range
armxstart=armx-range:armxend=armx+range
tblxstart=tblx-range:tblxend=tblx+range

OPEN "mk6results.txt" FOR OUTPUT AS #1
SCREEN 12:ProcTryLoops:PRINT " Done ":SLEEP
CLOSE #1
END'----------------------------------------------------------
SUB ProcTryLoops
    FOR a=astart TO aend STEP ink
        FOR b=bstart TO bend STEP ink
            FOR c=cstart TO cend STEP ink
                FOR d=dstart TO dend STEP ink
                    FOR e=estart TO eend STEP ink
                        FOR f=fstart TO fend STEP ink
                            FOR h=hstart TO hend STEP ink
                                FOR armx=armxstart TO armxend STEP ink
                                    FOR tblx=tblxstart TO tblxend STEP ink
                                        ProcVertexCoords'knt+=1:
                                    NEXT
                                NEXT
                            NEXT
                        NEXT
                    NEXT
                NEXT
            NEXT
            PRINT "b=";b;
        NEXT
        PRINT "a=";a
    NEXT
END SUB'--------------------------------------------------------
SUB ProcVertexCoords
    DIM AS DOUBLE dia1,dia2,dia3,dia4,diff
    '1st triangle and circumscribed circle diameter
    tbl=tblx+tinc1:arm=armx+ainc1:ProcForward:x1=x:y1=y'3 vertices
    tbl=tblx+tinc4:arm=armx+ainc4:ProcForward:x4=x:y4=y
    tbl=tblx+tinc7:arm=armx      :ProcForward:x7=x:y7=y
    s1=SQR(((x1-x4)^2)+((y1-y4)^2)):s2=SQR(((x4-x7)^2)+((y4-y7)^2)):s3=SQR(((x7-x1)^2)+((y7-y1)^2)) '3 sides
    dia1  =(2*s1*s2*s3)/SQR((s1+s2+s3)*(-s1+s2+s3)*(s1-s2+s3)*(s1+s2-s3))  'diameter
    'PRINT #1, 8888,a,b,dia1
    'PRINT #1, a,b,c,d,e,f,h,armx,tblx
    'PRINT #1, x1, y1, x4, y4, x7, y7
    
    IF dia1 > 155.67 AND dia1 < 155.69 THEN
        '2nd triangle and circumscribed circle diameter
        tbl=tblx+tinc2:arm=armx+ainc2:ProcForward:x2=x:y2=y
        tbl=tblx      :arm=armx+ainc5:ProcForward:x5=x:y5=y
        tbl=tblx+tinc8:arm=armx+ainc8:ProcForward:x8=x:y8=y
        s4=SQR(((x2-x5)^2)+((y2-y5)^2)):s5=SQR(((x5-x8)^2)+((y5-y8)^2)):s6=SQR(((x8-x2)^2)+((y8-y2)^2))
        dia2 =(2*s4*s5*s6)/SQR((s4+s5+s6)*(-s4+s5+s6)*(s4-s5+s6)*(s4+s5-s6))
    ELSE
        EXIT SUB
    END IF
    
    IF dia2 > 155.67 AND dia2 < 155.69 THEN
        '3rd triangle and circumscribed circle diameter
        tbl=tblx+tinc3:arm=armx+ainc3:ProcForward:x3=x:y3=y
        tbl=tblx+tinc6:arm=armx+ainc6:ProcForward:x6=x:y6=y
        tbl=tblx+tinc9:arm=armx+ainc9:ProcForward:x9=x:y9=y
        s7=SQR(((x3-x6)^2)+((y3-y6)^2)):s8=SQR(((x6-x9)^2)+((y6-y9)^2)):s9=SQR(((x9-x3)^2)+((y9-y3)^2))
        dia3 =(2*s7*s8*s9)/SQR((s7+s8+s9)*(-s7+s8+s9)*(s7-s8+s9)*(s7+s8-s9))
    ELSE
        EXIT SUB
    END IF
    
    IF dia3 > 155.67 AND dia3 < 155.69 THEN
        '4th triangle and circumscribed circle diameter
        tbl=tblx+tinc0:arm=armx+ainc0:ProcForward:x0=x:y0=y
        s10=SQR(((x0-x3)^2)+((y0-y3)^2)):s11=SQR(((x3-x7)^2)+((y3-y7)^2)):s12=SQR(((x7-x0)^2)+((y7-y0)^2))
        dia4 =(2*s10*s11*s12)/SQR((s10+s11+s12)*(-s10+s11+s12)*(s10-s11+s12)*(s10+s11-s12))
    ELSE 
        EXIT SUB
    END IF
    
    IF dia4 > 155.67 AND dia4 < 155.69 THEN
        ProcCircumscribed
    ELSE 
        EXIT SUB
    END IF
END SUB'--------------------------------------------------------
SUB ProcForward ' Screwlength to Vertex coords
    DIM AS DOUBLE g,aH,aX,aF,aG,aTBL,aARM'variables
    aB=ACOS(((a^2)+(c^2)-(b^2))/(2*a*c)):aC=ACOS(((a^2)+(b^2)-(c^2))/(2*a*b))
    aH=ACOS(((f*f)+(e*e)-(h*h))/(2*f*e))
    aTBL=ACOS(((d*d)+(b*b)-(tbl*tbl))/(2*d*b))     '<<<<<<< tbl
    aARM=ACOS(((c*c)+(e*e)-(arm*arm))/(2*c*e))     '<<<<<<< arm
    aG  =aB-aARM+aH
    aF  =ATN((f*(SIN(aG)))/(a-(f*(COS(aG)))))
    aX  =aTBL-(aC-aF)-0.4837
    g   =SQR((f*f)+(a*a)-(2*f*a*(COS (aG))))
    x   =(g*SIN(aX))+110                           'x >>>>>>>
    y   =572-(g*COS(aX))                           'y >>>>>>>
END SUB'---------------------------------------------------------

SUB ProcCircumscribed  'Center coordinates
    DIM AS DOUBLE Xu1,Yu1,d1,Xu2,Yu2,d2,Xu3,Yu3,d3,Xu4,Yu4,d4,xmax,xmin,ymax,ymin,xdiff,ydiff      
    
    d1=2*(x1*(y4-y7)+x4*(y7-y1)+x7*(y1-y4))'1 4 7                                ' <<<<< 10 Vertex coords for 4 triangles
    Xu1 =(((x1^2+y1^2)*(y4-y7))+((x4^2+y4^2)*(y7-y1))+((x7^2+y7^2)*(y1-y4)))/d1
    Yu1 =(((x1^2+y1^2)*(x7-x4))+((x4^2+y4^2)*(x1-x7))+((x7^2+y7^2)*(x4-x1)))/d1
    
    d2=2*(x2*(y5-y8)+x5*(y8-y2)+x8*(y2-y5))'2 5 8
    Xu2 =(((x2^2+y2^2)*(y5-y8))+((x5^2+y5^2)*(y8-y2))+((x8^2+y8^2)*(y2-y5)))/d2
    Yu2 =(((x2^2+y2^2)*(x8-x5))+((x5^2+y5^2)*(x2-x8))+((x8^2+y8^2)*(x5-x2)))/d2
    
    d3=2*(x3*(y6-y9)+x6*(y9-y3)+x9*(y3-y6))'3 6 9
    Xu3 =(((x3^2+y3^2)*(y6-y9))+((x6^2+y6^2)*(y9-y3))+((x9^2+y9^2)*(y3-y6)))/d3
    Yu3 =(((x3^2+y3^2)*(x9-x6))+((x6^2+y6^2)*(x3-x9))+((x9^2+y9^2)*(x6-x3)))/d3
    
    d4=2*(x0*(y3-y7)+x3*(y7-y0)+x7*(y0-y3))'0 3 7
    Xu4 =(((x0^2+y0^2)*(y3-y7))+((x3^2+y3^2)*(y7-y0))+((x7^2+y7^2)*(y0-y3)))/d4
    Yu4 =(((x0^2+y0^2)*(x7-x3))+((x3^2+y3^2)*(x0-x7))+((x7^2+y7^2)*(x3-x0)))/d4
    
    IF Xu1 > Xu2 THEN xmax=Xu1:xmin=Xu2 ELSE xmax=Xu2:xmin=Xu1
    IF Xu3 > xmax THEN xmax=Xu3
    IF Xu3 < xmin THEN xmin=Xu3
    IF Xu4 > xmax THEN xmax=Xu4
    IF Xu4 < xmin THEN xmin=Xu4
    xdiff=xmax-xmin
    IF Yu1 > Yu2 THEN ymax=Yu1:ymin=Yu2 ELSE ymax=Yu2:ymin=Yu1
    IF Yu3 > ymax THEN ymax=Yu3
    IF Yu3 < ymin THEN ymin=Yu3
    IF Yu4 > ymax THEN ymax=Yu4
    IF Yu4 < ymin THEN ymin=Yu4
    ydiff=ymax-ymin
    IF xdiff<0.12 AND ydiff<0.12 THEN PRINT #1, a,b,c,d,e,f,h,armx,tblx:PRINT #1, Xu1,Yu1:PRINT #1,Xu2,Yu2:PRINT #1,Xu3,Yu3:PRINT #1,Xu4,Yu4
END SUB'----------------------------------------------------------
