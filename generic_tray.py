#!/usr/bin/python3

helpTxt = """
ALL UNITS ARE MILLIMETERS

The tray will be specified by the dimensions of each row & column,
provided as two arrays.  If your bin is going to look like this:

               ----------------------------------
               |           |     |              |
        60mm   |           |     |              |
               ----------------------------------
               |           |     |              |
        60mm   |           |     |              |
               ----------------------------------
               |           |     |              |
               |           |     |              |
       100mm   |           |     |              |
               ----------------------------------
        30mm   |           |     |              |
               ----------------------------------
                   40mm      25mm     70mm


Use the following command-line call to generate the above (x-axis always first)
You can add extra -- args to adjust floor- and wall-thickness, as well as
bin depth and how rounded you want the bins to be on the bottom.
Make sure there are no spaces before or after the commas in the lists:

   python generic_tray.py [40,25,70] [30,100,60,60]
   python generic_tray.py [40,25,70] [30,100,60,60] --depth=55 --wall=1.5

Type the following to see all options:

   python generic_tray.py --help

*** IF YOU DON'T WANT/CAN'T DO THE COMMAND-LINE THING, then follow the
    directions towards where CLI_OPTIONS and CLI_ARGS are set.  You can
    remove all the CLI_* code and simply hardcode the tray size into this
    script so you only have to double-click it to generate the SCAD file.

"""
from solid import *
from solid.utils import *
from ast import literal_eval
from math import sqrt, pi
import os
import time
import sys
import optparse
import subprocess

try:
    # http://www.thingiverse.com/thing:230596
    # Use the printer calibration cubes to compute the 8 values and put
    # them into caldata.py.
    #
    # x/y/zScale is the raw multiplicative error in the printer in each dim
    # xi/yiAdj is inner (slot) c-values, xo/yoAdj is outer (peg) c-values
    # The *Adj values presumably have to do with filament width and shells
    from caldata import xScale,yScale,zScale, xiAdj,yiAdj,zAdj, xoAdj,yoAdj
except:
    print('*' * 80)
    print('***NO CALIBRATION DATA FOUND.  ASSUMING PERFECT CALIBRATION')
    print('*' * 80)
    xScale,yScale,zScale = 1.0, 1.0, 1.0
    xiAdj, yiAdj, zAdj = 0.0, 0.0, 0.0
    xoAdj, yoAdj = 0.0, 0.0

################################################################################
# This will create a plug that can be subtracted from the tray frame/box
def createSubtractSlot(offsetX, offsetY, sizeX, sizeY, binDepth,
                       roundDepth=15, trayFloor=1.5):

    sizeX = float(sizeX)
    sizeY = float(sizeY)

    # If round-depth is zero, it's just a square plug
    if roundDepth <= 0:
        return translate( [offsetX, offsetY, trayFloor]) \
            (
                cube( [sizeX, sizeY, depth * 1.1])
            )

    # Create 1:1 aspect, then stretch the whole thing at once
    # Prism sitting with corner at origin
    fullPrism = cube( [sizeX, sizeX, depth * 1.1])


    # Prism translated in the z-dir by the roundDepth
    partPrism = translate( [0, 0, roundDepth]) \
        (
            cube( [sizeX, sizeX, depth * 1.1])
        )

    # Start by creating a sphere in the center, scale it in y- an z-, then
    # translate it to the bottom of the partPrism
    sphereRad = sqrt(2) * sizeX / 2.0
    sphereScaleZ = roundDepth / sphereRad

    theSphere = translate( [sizeX / 2.0, sizeX / 2.0, roundDepth]) \
        (
            scale( [1, 1, sphereScaleZ]) \
            (
                sphere(sphereRad)
            )
        )

    return translate( [offsetX, offsetY, trayFloor]) \
        (
            scale( [1, sizeY / sizeX, 1]) \
            (
                intersection() \
                (
                    fullPrism,
                    union() \
                    (
                        partPrism,
                        theSphere
                    )
                )
            )
        )

################################################################################
def computeSlotVolume(xsz, ysz, depth, rdepth):
    """
    UPDATED:  Now we do this calculation exactly.  We compute the volume of the
    prism and then compute the volume of the hemisphere.  The hemisphere is
    complex because it's actually a hemisphere intersected with a square peg.

    We do the volume calculation by computing the volume of the full hemisphere
    and then subtracting the volume of the four "spherical caps".   At once we
    have that, we scale the volume by both the y-scale and z-scale.

    From http://en.wikipedia.org/wiki/Spherical_cap the volume of a spherical
    cap is:

       pi * h * (3a*a + h*h) / 6

    "h" is the height of the cap which is the radius of sphere minus x/2
    "a" is the radius of the base of the cap, which is just x/2

    Don't forget to cut the resultant sph cap volume in half, because we're
    only removing half of a cap (because it's from half a hemisphere)
    """

    #
    sphereRad = sqrt(2) * xsz / 2.0
    sphereScaleY = ysz / xsz
    sphereScaleZ = rdepth / sphereRad

    a = xsz / 2.0
    h = sphereRad - a
    oneCapFullVol = pi * h * (3 * a * a + h * h) / 6.0

    fullSphereVol = 4.0 * pi * sphereRad ** 3 / 3.0

    roundVol_mm3 = (fullSphereVol - 4 * oneCapFullVol) / 2.0
    roundVol_mm3 *= sphereScaleY
    roundVol_mm3 *= sphereScaleZ

    prismTop_mm3 = (depth - rdepth) * xsz * ysz
    totalVol_mm3 = prismTop_mm3 + roundVol_mm3

    # Now convert to both cups and mL (imperial and metric)
    totalVol_cm3 = totalVol_mm3 / 10 ** 3
    totalVol_mL = totalVol_cm3           # 1 cm3 == 1 mL  !
    totalVol_cups = totalVol_mm3 / 236588.

    return [totalVol_mL, totalVol_cups]

def createTray(xlist, ylist, dep, rdep=15, wall=1.5, floor=1.5):

    # Create all the slots to be subtracted from the frame of the tray.
    slots = []
    xOff = wall
    yOff = wall

    for ysz in ylist:
        xOff = wall
        for xsz in xlist:
            slots.append(createSubtractSlot(xOff, yOff, xsz, ysz, dep, rdep, floor))
            xOff += wall + xsz
        yOff += wall + ysz

    # The loops leave xOff & yOff at the upper-left corner of the tray.  Perfect!
    totalWidth = xOff
    totalHeight = yOff

    # We have a list of "slots", unionize them.
    allStuffToSubtract = union() (*slots)

    # Create the prism from which the slots will be subtracted
    trayBasePrism = cube( [totalWidth, totalHeight, floor + depth])

    # Finally, create the object and scale by the printer-calibration data
    return [totalWidth, totalHeight,
            scale( [xScale, yScale, zScale]) \
                (
                    difference() \
                    (
                        trayBasePrism,
                        allStuffToSubtract
                    )
                )]

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="%prog [options]\n")
    parser.add_option("--depth", dest="depth", default=38.0, type="float",
                      help="Depth of the tray above floor (mm, default 40)")
    parser.add_option("--floor", dest="floor", default=1.0, type="float",
                      help="Thickness of tray floor (mm, default 1.0)")
    parser.add_option("--wall", dest="wall", default=1.5, type="float",
                      help="Thickness of walls (mm, default 0.5)")
    parser.add_option("--round", dest="rdepth", default=30, type="float",
                      help="Height of tapered bottom (mm, default 15)")
    parser.add_option("--outfile", dest="outfile", default='', type="str",
                      help="The output name of the resultant file")
    (CLI_OPTIONS, CLI_ARGS) = parser.parse_args()


    if len(CLI_ARGS) < 2:
        print('***ERROR: Must provide the bin sizes as two arrays')
        print(helpTxt)
        exit(0)

    # If you want to specify the args directly here instead of using the CLI,
    # comment/delete the above CLI_ARGS clause and then set the following
    # variables directly
    floor = CLI_OPTIONS.floor
    wall = CLI_OPTIONS.wall
    depth = CLI_OPTIONS.depth
    rdepth = CLI_OPTIONS.rdepth
    fname = CLI_OPTIONS.outfile
    xsizes = [float(x) for x in literal_eval(CLI_ARGS [0])] # should be a list
    ysizes = [float(x) for x in literal_eval(CLI_ARGS [1])] # should be a list

    # Example for replacing the above lines if you don't want to use CLI
    #floor  = 1.5
    #wall   = 1.5
    #depth  = 40
    #rdepth = 15
    #xsizes = [30,45,60]
    #ysizes = [50,50,50,50]

    if rdepth > depth - 5:
        print('***Warning:  round depth needs to be smaller or equal to bin depth')
        ok = input('Shorten round depth? [Y/n]: ')
        if ok.lower().startswith('n'):
            print('Aborting...')
            exit(1)
        else:
            rdepth = max(depth - 5, 0)

    xszStrs = [str(int(x)) for x in xsizes]
    yszStrs = [str(int(y)) for y in ysizes]

    print('Floor:  ', floor, 'mm')
    print('Wall:   ', wall, 'mm')
    print('Depth:  ', depth, 'mm')
    print('Round:  ', rdepth, 'mm')
    print('Widths: ', xsizes, 'mm')
    print('Heights:', ysizes, 'mm')

    # If you don't override fname, the scad file will automatically be named
    if not fname:
        fname = 'tray_%s_by_%s.scad' % ('x'.join(xszStrs), 'x'.join(yszStrs))

    print('Writing to OpenSCAD file:', fname)

    # Now tell solid python to create the .scad file
    twid,thgt,trayObj = createTray(xsizes, ysizes, depth, rdepth, wall, floor)
    scad_render_to_file(trayObj, fname, file_header='$fn=64;')

    ################################################################################
    # EVERYTHING BELOW THIS LINE IS SIMPLY FOR PRINTING ASCII DIAGRAMS
    # Print some useful info
    # xOff and yOff end in the upper-right corner... which is the tray size!
    print('Tray size is: %0.2fmm by %0.2fmm' % (twid, thgt))
    print('              (%0.2fcm by  %0.2fcm)' % (twid / 10, thgt / 10))

    # The diagram will be approximately 72 chars wide by 48 chars tall
    # Console letters are about 1.5 times taller than they are wide
    totalCharsWide = 82
    totalCharsHigh = float(thgt) / float(twid) * (totalCharsWide / 2.0)

    def getCharsWide(mm):
        maxInternalChars = totalCharsWide - (len(xsizes) + 1)
        return max(10, int(maxInternalChars * float(mm) / float(twid)))

    def getCharsHigh(mm):
        maxInternalChars = totalCharsHigh - (len(ysizes) + 1)
        return max(2, int(maxInternalChars * float(mm) / float(thgt)))

    xchars = [getCharsWide(x) for x in xsizes]
    ychars = [getCharsHigh(y) for y in ysizes]

    print('')
    print('')

    wCharsTotal = sum(xchars) + len(xchars) + 1
    hCharsTotal = sum(ychars) + len(ychars) + 1

    vertLine = ' ' * 10 + '-' * wCharsTotal + '\n'
    sys.stdout.write(vertLine)
    for j in range(len(ysizes)):
        # Acually do the y-values in reverse since printing to console happens
        # in the negative y-direction.
        revj = len(ysizes) - j - 1

        for jc in range(ychars [revj]):
            yhgt = ychars [revj]
            if jc == yhgt / 2:
                sys.stdout.write( ('%0.1f mm' % ysizes [revj]).center(10) + '|')
            else:
                sys.stdout.write(' ' * 10 + '|')

            for i in range(len(xsizes)):
                mL,cups = computeSlotVolume(xsizes [i], ysizes [revj], depth, rdepth)

                if jc == yhgt / 2 - 1:
                    sys.stdout.write( ('%0.2f cups' % cups).center(xchars [i]))
                elif jc == yhgt / 2:
                    sys.stdout.write( ('%0.1f mL' % mL).center(xchars [i]))
                else:
                    sys.stdout.write(' ' * (xchars [i]))
                sys.stdout.write('|')
            sys.stdout.write('\n')

        sys.stdout.write(vertLine)

    sys.stdout.write('\n')
    sys.stdout.write(' ' * 10)
    for i in range(len(xsizes)):
        sizeStr = '%0.1f mm' % xsizes [i]
        sys.stdout.write(sizeStr.center(xchars [i] + 1))
    sys.stdout.write('\n\n')

    print('Total Width  (with walls):  %3.1f mm   (%3.2f cm)' % (twid,  twid / 10.0))
    print('Total Height (with walls):  %3.1f mm   (%3.2f cm)' % (thgt,  thgt / 10.0))
    print('Total Depth  (with floor):  %3.1f mm   (%3.2f cm)' % (depth + floor, (depth + floor) / 10.0))
    print('')

    ok = input('Generate STL file? (this can take a few minutes) [N/y]: ')
    if ok.lower().startswith('y'):
        stlname = fname + '.stl'
        print('Converting to STL file:', stlname,)
        proc = subprocess.Popen('openscad -o "%s" "%s"' % (stlname, fname), shell=True)
        while proc.poll() == None:
            time.sleep(0.1)
