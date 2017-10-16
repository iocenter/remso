WELSPECS
--     WELL   GROUP LOCATION  BHP   PI
--     NAME   NAME    I  J   DEPTH DEFN
       'i1'   'grp1'     7  7    1*  'WAT'  /
       'i2'   'grp1'     20  1    1*  'WAT'  /
       'p1'   'grp1'     1  1    1*  'WAT'  /
       'p2'   'grp1'     1  10    1*  'WAT'  /
       'p3'   'grp1'     20  10    1*  'OIL'  /
       'p4'   'grp1'     14  4    1*  'WAT'  /
       'p5'   'grp1'     8  2    1*  'OIL'  /

 / 

COMPDAT
--     WELL     -LOCATION- OPEN/ SAT CONN  WELL
--     NAME     I  J K1 K2 SHUT  TAB FACT  DIAM
       'i1'     7  7  1  5  'OPEN' 2*       0.3/
       'i2'     20  1  1  5  'OPEN' 2*       0.3/
       'p1'     1  1  1  5  'OPEN' 2*       0.3/
       'p2'     1  10  1  5  'OPEN' 2*       0.3/
       'p3'     20  10  1  5  'OPEN' 2*       0.3/
       'p4'     14  4  1  5  'OPEN' 2*       0.3/
       'p5'     8  2  1  5  'OPEN' 2*       0.3/

 / 

-- CONTROL INTERVAL 1
WCONPROD
--      WELL     OPEN/  CNTL   OIL  WATER   GAS  LIQU   RES   BHP
--      NAME     SHUT   MODE  RATE   RATE  RATE  RATE  RATE
        'p1'    'OPEN'  'BHP'   5*                            180.08/
        'p2'    'OPEN'  'BHP'   5*                            180.019/
        'p3'    'OPEN'  'ORAT'  5  5*    /
        'p4'    'OPEN'  'BHP'   5*                            180.256/
        'p5'    'OPEN'  'ORAT'  5  5*    /
/ 

WCONINJE
--      WELL     INJ   OPEN/   CNTL    FLOW   RES    BHP
--      NAME    TYPE   SHUT    MODE    RATE   RATE  TARGET
        'i1'   'WATER' 'OPEN' 'RATE'   300/
        'i2'   'WATER' 'OPEN' 'RATE'   150/
/

TSTEP
30 2 4 8 16 32 50 50 50 50 50 50 50 50 50 34.5 
/

-- CONTROL INTERVAL 2
WCONPROD
--      WELL     OPEN/  CNTL   OIL  WATER   GAS  LIQU   RES   BHP
--      NAME     SHUT   MODE  RATE   RATE  RATE  RATE  RATE
        'p1'    'OPEN'  'BHP'   5*                            182.332/
        'p2'    'OPEN'  'BHP'   5*                            181.283/
        'p3'    'OPEN'  'ORAT'  5  5*    /
        'p4'    'OPEN'  'BHP'   5*                            183.6/
        'p5'    'OPEN'  'ORAT'  5  5*    /
/ 

WCONINJE
--      WELL     INJ   OPEN/   CNTL    FLOW   RES    BHP
--      NAME    TYPE   SHUT    MODE    RATE   RATE  TARGET
        'i1'   'WATER' 'OPEN' 'RATE'   200/
        'i2'   'WATER' 'OPEN' 'RATE'   150/
/

TSTEP
1 2 4 8 16 32 50 50 50 50 50 50 50 50 50 34.5 
/

-- CONTROL INTERVAL 3
WCONPROD
--      WELL     OPEN/  CNTL   OIL  WATER   GAS  LIQU   RES   BHP
--      NAME     SHUT   MODE  RATE   RATE  RATE  RATE  RATE
        'p1'    'OPEN'  'BHP'   5*                            185.718/
        'p2'    'OPEN'  'BHP'   5*                            189.806/
        'p3'    'OPEN'  'ORAT'  5  5*    /
        'p4'    'OPEN'  'BHP'   5*                            185.611/
        'p5'    'OPEN'  'ORAT'  5  5*    /
/ 

WCONINJE
--      WELL     INJ   OPEN/   CNTL    FLOW   RES    BHP
--      NAME    TYPE   SHUT    MODE    RATE   RATE  TARGET
        'i1'   'WATER' 'OPEN' 'RATE'   200/
        'i2'   'WATER' 'OPEN' 'RATE'   150/
/

TSTEP
1 2 4 8 16 32 50 50 50 50 50 50 50 50 50 34.5 
/

-- CONTROL INTERVAL 4
WCONPROD
--      WELL     OPEN/  CNTL   OIL  WATER   GAS  LIQU   RES   BHP
--      NAME     SHUT   MODE  RATE   RATE  RATE  RATE  RATE
        'p1'    'OPEN'  'BHP'   5*                            186.489/
        'p2'    'OPEN'  'BHP'   5*                            180/
        'p3'    'OPEN'  'ORAT'  5  5*    /
        'p4'    'OPEN'  'BHP'   5*                            182.09/
        'p5'    'OPEN'  'ORAT'  5  5*    /
/ 

WCONINJE
--      WELL     INJ   OPEN/   CNTL    FLOW   RES    BHP
--      NAME    TYPE   SHUT    MODE    RATE   RATE  TARGET
        'i1'   'WATER' 'OPEN' 'RATE'   200/
        'i2'   'WATER' 'OPEN' 'RATE'   150/
/

TSTEP
1 2 4 8 16 32 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 49.5 
/

-- CONTROL INTERVAL 5
WCONPROD
--      WELL     OPEN/  CNTL   OIL  WATER   GAS  LIQU   RES   BHP
--      NAME     SHUT   MODE  RATE   RATE  RATE  RATE  RATE
        'p1'    'OPEN'  'BHP'   5*                            186.352/
        'p2'    'OPEN'  'BHP'   5*                            180/
        'p3'    'OPEN'  'ORAT'  5  5*    /
        'p4'    'OPEN'  'BHP'   5*                            182.53/
        'p5'    'OPEN'  'ORAT'  5  5*    /
/ 

WCONINJE
--      WELL     INJ   OPEN/   CNTL    FLOW   RES    BHP
--      NAME    TYPE   SHUT    MODE    RATE   RATE  TARGET
        'i1'   'WATER' 'OPEN' 'RATE'   200/
        'i2'   'WATER' 'OPEN' 'RATE'   150/
/

TSTEP
1 2 4 8 16 32 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 32 
/

