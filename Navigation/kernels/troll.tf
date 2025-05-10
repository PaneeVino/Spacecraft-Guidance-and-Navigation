KPL/FK
 
   FILE: troll.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2024-11-20T12:37:44
   PINPOINT DEFINITIONS FILE: troll.def
   PINPOINT PCK FILE:         pck00010.tpc
   PINPOINT SPK FILE:         troll.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'TROLL'
   NAIF_BODY_CODE                      += 399101
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame TROLL_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame TROLL_TOPO is centered at the
      site TROLL, which has Cartesian coordinates
 
         X (km):                  0.1974141000000E+04
         Y (km):                  0.8743900000000E+02
         Z (km):                 -0.6045333000000E+04
 
      and planetodetic coordinates
 
         Longitude (deg):         2.5360971631568
         Latitude  (deg):       -72.0119786290199
         Altitude   (km):         0.1297635025147E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame IAU_EARTH.
 
 
\begindata
 
   FRAME_TROLL_TOPO                    =  1399101
   FRAME_1399101_NAME                  =  'TROLL_TOPO'
   FRAME_1399101_CLASS                 =  4
   FRAME_1399101_CLASS_ID              =  1399101
   FRAME_1399101_CENTER                =  399101
 
   OBJECT_399101_FRAME                 =  'TROLL_TOPO'
 
   TKFRAME_1399101_RELATIVE            =  'IAU_EARTH'
   TKFRAME_1399101_SPEC                =  'ANGLES'
   TKFRAME_1399101_UNITS               =  'DEGREES'
   TKFRAME_1399101_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399101_ANGLES              =  (   -2.5360971631568,
                                            -162.0119786290199,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file troll.def
--------------------------------------------------------------------------------
 
begindata
 
      SITES                      += 'TROLL'
      TROLL_FRAME                  = 'IAU_EARTH'
      TROLL_CENTER                 = 399
      TROLL_IDCODE                 = 399101
      TROLL_BOUNDS                 = ( @1950-JAN-01/00:00,  @2050-JAN-01/00:00 )
      TROLL_XYZ                    = (  1974.141
                                       87.439
                                       -6045.333
                                   )
      TROLL_UP      = 'Z'
      TROLL_NORTH   = 'X'
 
begintext
 
begintext
 
[End of definitions file]
 
