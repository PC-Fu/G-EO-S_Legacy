#!python
# Borehole meshing
# Lee Aarons & Ryan Vignes
# 

# Can be run as script or interactively (within Cubit GUI) -- make designated changes below regarding input file (at the start of the "Read input file" section)

# Creates a rectangular prism cut through the center by a borehole (oriented in the z-direction)
# Ellipsoid perforations are carved out of the side of the borehole, at locations dictated by phase angle and spacing
# An ellipsoid shell surrounding each perforation is cut out and reattached. These shells allow for the mesh size near and far from the perforation to be different and generally reduces (and often eliminates) meshing errors and bad elements.

# --------------
# Initialize
# --------------

import sys
sys.path.append("/usr/gapps/cubit/linux64.14.0/bin")
import cubit

cubit.init([""])

c = cubit.cmd
cc = cubit
c('undo off') # turn undo off to speed up, may want to turn on in testing, etc.

# -----------------------------
# Read input file
# -----------------------------

infname = sys.argv[1]                  # input file name, comment out if running interactively and use next line
# infname = "borehole75.txt"             # input file name if running interactively
inputfile = open(infname,"r")
outname = inputfile.readline()         # abaqus file to output (generally .inp file)
cubname = inputfile.readline()         # cubit save file (generally .cub file)
errname = inputfile.readline()         # error logging file (generally .err file)
c('Logging Errors On File "%s"' % (errname))
xsize = float(inputfile.readline())    # full system size (x-dir)
ysize = float(inputfile.readline())    # full system size (y-dir)
zsize = float(inputfile.readline())    # full system size (z-dir)
holerad = float(inputfile.readline())  # borehole radius
perfwid = float(inputfile.readline())  # perforation width (diameter) at entrance
perflen = float(inputfile.readline())  # perforation length
phase = float(inputfile.readline())    # perforation phasing (angle in degrees)
spacing = float(inputfile.readline())  # perforation spacing
numperfs = int(inputfile.readline())   # number of perforations
firstperf = float(inputfile.readline()) # depth within borehole where perfs start
coarsesize = float(inputfile.readline())    # mesh size (avg. hex vol.) "far" from perfs 
finesize = float(inputfile.readline())  # mesh size (avg. hex vol.) around perfs
uniteatend = int(inputfile.readline())  # unite volumes into 1? -- 1 if yes, no otherwise
addhalf = 2
addhalf = int(inputfile.readline())    # use webcuts at half-angles -- 1 == yes, 2 == auto, no otherwise
inputfile.close()

# Over the course of generating the geometry, it is cut into many pieces. Set uniteatend to 1 to unite
# the pieces into a single volume after meshing. Advantages of doing so are (1) Simpler/smaller abaqus 
# file (assuming element blocks are kept), (2) Easier to visualize, and (3) Should hexes with negative 
# Jacobians still exist, smoothing after uniting sometimes helps (after latest changes to script, this may
# be unnecessary). Disadvantages of uniting are (1) It takes a long time and (2) Errors or crashing may 
# occur. The current recommendation is to set uniteatend to 0.

# One set of cutting is performed along the length of borehole at angles equal to the perforation angles. 
# If the phase angle is p degrees and there are n perforations, n such cuts are made separated by p 
# degrees. For example, if the phase angle is 60 degrees, assuming there are at least 6 perforations, cuts
# will be performed at 0, 60, 120, 180, 240, and 300 degrees. If addhalf is set to 1, cutting is also 
# performed at the half-angle (every p/2 degrees). For situations in which phase angle = 360/m, with m =
# integer, addhalf should be set to 1. Otherwise, if there are enough perfs to cover 360 degrees addhalf 
# should (probably) be equal to 0. Setting addhalf to 2 will let the script choose whether to use 
# half-angle cuts automaticlly -- this is default and recommended.

ellipse_xcoords   = [0.0, 0.0, perfwid/2.]  # List of x-coordinates forming perf, in ellipse's reference frame
ellipse_ycoords   = [0.0, perflen, 0.0]  # List of y-coordinates forming perf, in ellipse's reference frame
ellipse_zcoords   = [0.0, 0.0, 0.0]  # List of z-coordinates forming perf, in ellipse's reference frame 

xz_hole_fact      = 2.5   # Scale factor on ellipse region surrounding perfs for local mesh refinement
y_hole_fact       = 1.5   # Scale factor on ellipse region surrounding perfs for local mesh refinement

# -----------------------------
# Turn off visualaization (faster run time) -- uncomment if running interactively and want speed up
# -----------------------------
# c('geometry visibility off')
# c('mesh visibility off')
# c('undo off')

# -----------------------------
# Subroutines for mesher
# -----------------------------
def webcut_group_with_plane(grp_name,plane_dir,plane_offset):
   c('${temp_vol_id = Id("volume")}')
   c('webcut volume in %s with plane xplane rotate %s about z center 0 0 0 noimprint nomerge' % (grp_name, plane_dir) )
   c('group "%s" add volume {temp_vol_id+1}' % (grp_name) )
   c('webcut volume in %s with plane zplane offset %f' % (grp_name,plane_offset) )
   c('group "%s" add volume {temp_vol_id+2}  {temp_vol_id+3}' % grp_name )

# -----------------------------
# Create geometry
# -----------------------------
# Main bore hole, axial aligned with z-axis

c('brick x %s y %s z %s' % (xsize, ysize, zsize))
fullsystem = cc.get_last_id("volume")

# ---- Create perf geometry, will be replicated for all perfs
vrtcs = []
for jj in range( len(ellipse_xcoords) ):
   c('create vertex x %s y %s z %s' % (ellipse_xcoords[jj], ellipse_ycoords[jj], ellipse_zcoords[jj]) )
   vrtcs.append(cc.get_last_id("vertex"))

c('create curve vertex %s vertex %s vertex %s ellipse first angle 0 last angle 170' % (vrtcs[1], vrtcs[2], vrtcs[0]))
ecurve = cc.get_last_id("curve")
c('delete vertex all')
c('sweep curve %s yaxis angle 360 make_solid' % (ecurve))

hole_cutout = cc.get_last_id("volume")

# ---- Scale hole shape to create a slightly larger version for
# ---- mesh refinement surrounding the hole

c('volume %s scale x %s y %s z %s copy' % (hole_cutout, xz_hole_fact, y_hole_fact, xz_hole_fact) )
scaled_hole = cc.get_last_id("volume")

c('move volume %s y %s' % (hole_cutout, holerad))
c('move volume %s y %s' % (scaled_hole, holerad))

# create borehole

holelen = zsize
c('create cylinder height %s radius %s' % (holelen, holerad)) # (bore_height, bore_in_rad))
cylvol = cc.get_last_id("volume")
c('move volume %s z %s include_merged ' % (cylvol, (holelen - zsize)/2.)) # ((bore_height)/2) )
c('subtract volume %s from volume 1 %s %s' %(cylvol, hole_cutout, scaled_hole) )

# ---- Move hole geometry to new locations

hole_groups = []
zplace = firstperf - zsize/2.
angle = 360.0
cut_planes = []
cut_angles = []
hole_grpnums = []
bslc_grpnums = []
for jj in range(numperfs): 
   #  copy surrounding volume and remove it to create surface, then recopy the volume
   c('volume %s copy move x 0.0 y 0.0 z %s' % (scaled_hole, zplace))
   latest_vol = cc.get_last_id("volume")
   c('rotate volume %s angle %s about z include_merged' % (latest_vol, angle))
   c('subtract volume %s from %s imprint' % (latest_vol, fullsystem))
   c('volume %s copy move x 0.0 y 0.0 z %s' % (scaled_hole, zplace))
   latest_vol = cc.get_last_id("volume")
   c('rotate volume %s angle %s about z include_merged' % (latest_vol, angle))
   c('volume %s copy move x 0.0 y 0.0 z %s' % (hole_cutout, zplace))
   last_hole = cc.get_last_id("volume")
   c('rotate volume %s angle %s about z include_merged' % (last_hole, angle))
   c('subtract volume %s from volume %s' % (last_hole, latest_vol))
   cut_planes.append(zplace)
   cut_angles.append(angle)
   aa = "hole" + str(latest_vol)
   hole_groups.append(aa)
   c('group "%s" add volume %s' % (hole_groups[jj], latest_vol) )
   hole_grpnums.append(cc.get_last_id("group"))
   webcut_group_with_plane(hole_groups[jj], angle, zplace)
   zplace = zplace + spacing
   angle = angle + phase

c('delete volume %s' % (hole_cutout))
c('delete volume %s' % (scaled_hole) )

# ---- Divide bore hole into planes centered along hole height location

plane_grp = []
for jj in range( len(cut_planes) ):
   plane_grp.append("bore_plane" + str(jj))
   c('webcut volume %s with plane zplane offset %s noimprint nomerge' % (fullsystem, cut_planes[jj])) #main_bore_cyl, cut_planes[jj]) )
   bore_slice = cc.get_last_id("volume")
   c('group "%s" add volume %s' % (plane_grp[jj], bore_slice) )
   bslc_grpnums.append(cc.get_last_id("group"))

c('group "%s" add volume %s' % ("bore_plane" + str(jj+1), fullsystem))  #main_bore_cyl) )
plane_grp.append("bore_plane" + str(jj+1))
bslc_grpnums.append(cc.get_last_id("group"))

# ---- Slice bore planes along direction of holes

old_cut_angles = cut_angles
cut_angles = []
if addhalf == 2:
   if 360.%phase == 0.:
      addhalf = 1
   elif numperfs*phase < 360.:
      addhalf = 1
   else:
      addhalf = 0

for i in range(numperfs):
   angle = (i*phase) % 360.0 + 360.0
   cut_angles.append(angle)
   if addhalf == 1:
      angle = ((i+0.5)*phase) % 360.0 + 360.0
      cut_angles.append(angle)

angle_planes = list( set(cut_angles)) #side_hole_angle) )
for ii in range( len(plane_grp) ):
   for jj in range( len(angle_planes) ):
      c('webcut volume in %s with plane xplane rotate %s about z center 0 0 0 noimprint nomerge' % (plane_grp[ii], angle_planes[jj]) )
      current_vol = cc.get_last_id("volume")
      if(jj == 0):
         c('group "%s" add volume %s' % (plane_grp[ii], current_vol) ) 
      else:
         c('group "%s" add volume %s %s' % (plane_grp[ii], current_vol-1, current_vol) )

# -----------------------------
# Mesh geometry
# -----------------------------
c('imprint volume all')
c('merge volume all')

if uniteatend != 1:   # if not uniting volumes into one body, easier to assign surfaces to nodesets before meshing
   c('undo on')
   c('create cylinder height %s radius %s' % (holelen, holerad))
   newcyl = cc.get_last_id("volume")
   c('imprint volume all')
   c('merge volume all')
   boresurfaces = cc.get_relatives("volume",newcyl,"surface")
   c('undo')
   c('undo')
   c('nodeset 7 surface all in volume %s' % (newcyl))
   c('delete volume %s' % (newcyl))
   holesurfs = []
   for ii in hole_grpnums:
      vols = cc.get_group_volumes(ii)
      for vol in vols:
         surfs = cc.get_relatives("volume",vol,"surface")
         for surf in surfs:
            holesurfs.append(surf)
   surroundsurfs = []
   for ii in boresurfaces:
      surroundsurfs.append(ii)
      c('nodeset 7 surface %s' %(ii))
   for ii in bslc_grpnums:
      vols = cc.get_group_volumes(ii)
      for vol in vols:
         surfs = cc.get_relatives("volume",vol,"surface")
         for surf in surfs:
            surroundsurfs.append(surf)            
   for jj in holesurfs:
      if jj not in surroundsurfs:
         c('nodeset 8 surface %s' % (jj))
   c('nodeset 1 surface all with x_coord = %s' % (-xsize/2)) 
   c('nodeset 2 surface all with y_coord = %s' % (-ysize/2))
   c('nodeset 3 surface all with z_coord = %s' % (-zsize/2))
   c('nodeset 4 surface all with x_coord = %s' % (xsize/2)) 
   c('nodeset 5 surface all with y_coord = %s' % (ysize/2))
   c('nodeset 6 surface all with z_coord = %s' % (zsize/2))
   c('set print quality off')
   c('set continue meshing off')

c('set file overwrite on')
c('save as "' + cubname + '" Overwrite')  # save .cub file before meshing

for hole in hole_groups:
   c('volume in group %s size %s' % (hole, finesize) )

for bslice in plane_grp:
   c('volume in group %s size %s' % (bslice, coarsesize))

c('volume all scheme polyhedron')

for bslice in plane_grp:
   c('mesh volume in group %s' % (bslice))
   c('volume in group %s smooth scheme untangle beta 0.02 cpu 10' % (bslice) )
   c('smooth volume in group %s' % (bslice) )
#
for hole in hole_groups:
   c('mesh volume in group %s' % (hole))
   c('volume in group %s smooth scheme untangle beta 0.02 cpu 10' % (hole) )
   c('smooth volume in group %s' % (hole) )

if uniteatend == 1:  # if unite to make single volume, unite and assign nodesets
   c('unite all include_mesh')
   c('compress')
   c('volume 1 smooth scheme untangle beta 0.02 cpu 10')
   c('smooth volume 1')
   nsurfs = numperfs + 7
   for surfs in range(nsurfs):
      stype = cc.get_surface_type(surfs+1)
      if stype == 'cone surface':
         c('nodeset 7 surface %s' % (surfs+1))
      elif stype == 'spline surface':
         c('nodeset 8 surface %s' % (surfs+1))
      elif stype == 'plane surface': 
         snorm = cc.get_surface_normal(surfs+1)
         if (snorm[0] < -0.9):
            c('nodeset 1 surface %s' % (surfs+1))
         elif (snorm[0] > 0.9):
            c('nodeset 4 surface %s' % (surfs+1))
         elif (snorm[1] < -0.9):
            c('nodeset 2 surface %s' % (surfs+1))
         elif (snorm[1] > 0.9):
            c('nodeset 5 surface %s' % (surfs+1))
         elif (snorm[2] < -0.9):
            c('nodeset 3 surface %s' % (surfs+1))
         elif (snorm[2] > 0.9):
            c('nodeset 6 surface %s' % (surfs+1))
   c('nodeset 7 name "borehole"')
   c('nodeset 2 name "yneg"')
   c('nodeset 3 name "zneg"')
   c('nodeset 1 name "xneg"')
   c('nodeset 5 name "ypos"')
   c('nodeset 6 name "zpos"')
   c('nodeset 4 name "xpos"')
   c('nodeset 8 name "cav"')

mcmnd = 'export abaqus "' + outname + '" dimension 3 overwrite cubitids'
c(mcmnd)

cc.destroy()
