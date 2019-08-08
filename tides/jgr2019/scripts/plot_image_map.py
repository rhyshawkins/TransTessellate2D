
import sys
import numpy

import argparse

import matplotlib.pyplot as P
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable

def load_regression_points(filename):
    f = open(filename, 'r')
    lines = f.readlines()[1:]
    f.close()

    return zip(*map(lambda x: map(float, x.split())[:3], lines))

def load_combined_points(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    tide = []  #1
    land = []  #2
    sea = []   #3
    
    
    for line in lines[1:]:
        t = line.split()
        lon = float(t[0])
        lat = float(t[1])
        obstype = int(t[2])
        value = float(t[3])
        err = float(t[4])

        if obstype == 1:
            tide.append((lon, lat, value, err))
        elif obstype == 2:
            land.append((lon, lat, value, err))
        elif obstype == 3:
            sea.append((lon, lat, value, err))
        else:
            raise Exception('Unknown observation type %d' % obstype)


    return numpy.array(tide), numpy.array(land), numpy.array(sea)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')
    parser.add_argument('--subtract', type = str, default = None, help = 'Subtract from image')
    
    parser.add_argument('--flip', action = 'store_true', default = False, help = 'Vertically flip image')
    
    parser.add_argument('-x', '--minlon', type = float, default = -1.0, help = 'Min Lon')
    parser.add_argument('-X', '--maxlon', type = float, default = 1.0, help = 'Max Lon')
    parser.add_argument('-y', '--minlat', type = float, default = -1.0, help = 'Min Lat')
    parser.add_argument('-Y', '--maxlat', type = float, default = 1.0, help = 'Max Lat')

    parser.add_argument('--vmin', type = float, default = 0.0, help = 'V Min')
    parser.add_argument('--vmax', type = float, default = 1.0, help = 'V Max ')

    parser.add_argument('--colorbar', action = 'store_true', default = False, help = 'Show colour bar')
    parser.add_argument('--horizontal', action = 'store_true', default = False, help = 'Horizontal color bar')
    parser.add_argument('--clabel', type = str, default = None, help = 'Colour bar label')

    parser.add_argument('--pdf', type = str, default = None, help = 'PDF output')
    parser.add_argument('--png', type = str, default = None, help = 'PNG output')

    parser.add_argument('--noshow', action = 'store_true', default = False, help = 'No output')

    parser.add_argument('--invertcoastlines', action = 'store_true', default = False, help = 'Border colour')

    parser.add_argument('--parallels', type = int, default = 10, help = 'Parallel spacing')
    parser.add_argument('--meridians', type = int, default = 10, help = 'Meridians spacing')

    parser.add_argument('--cmap', type = str, default = 'RdBu', help = 'Colour map')

    parser.add_argument('--resolution', type = str, default = 'c', help = 'Resolution of boundaries')
    parser.add_argument('--fine', action = 'store_const', dest = 'resolution', const = 'f', help = 'Fine boundaries')
    parser.add_argument('--intermediate', action = 'store_const', dest = 'resolution', const = 'i', help = 'Intermediate resolution boundaries')


    parser.add_argument('--points', type = str, default = None, help = 'Show points')
    parser.add_argument('--point-filter', type = int, default = 0, help = 'Show only selected points')
    parser.add_argument('--tide', action = 'store_const', dest = 'point_filter', const = 1, help = 'Show  only tide points')
    parser.add_argument('--land', action = 'store_const', dest = 'point_filter', const = 2, help = 'Show  only land points')
    parser.add_argument('--sea', action = 'store_const', dest = 'point_filter', const = 3, help = 'Show  only sea points')

    
    args = parser.parse_args()
    
    image = numpy.loadtxt(args.input)

    rows, cols = image.shape

    if not (args.subtract is None):
        subtract = numpy.loadtxt(args.subtract)

        image = image - subtract
        

    if args.flip:
        image = numpy.flipud(image)

    print(numpy.min(image), numpy.max(image))

    dlon = (args.maxlon - args.minlon)/float(cols)
    dlat = (args.maxlat - args.minlat)/float(rows)
    
    lons = numpy.linspace(args.minlon, args.maxlon, cols + 1)[:-1] + dlon/2.0
    lats = numpy.linspace(args.minlat, args.maxlat, rows + 1)[:-1] + dlat/2.0

    fig, ax = P.subplots()
    #fig.set_tight_layout(True)

    m = Basemap(resolution=args.resolution,
                projection='merc',
                llcrnrlon = args.minlon,
                urcrnrlon = args.maxlon,
                llcrnrlat = args.minlat,
                urcrnrlat = args.maxlat,
                ax = ax)
    
    x, y = numpy.meshgrid(lons, lats)

    p = m.pcolor(x, y, image,
                 latlon = True,
                 cmap = args.cmap,
                 vmin = args.vmin,
                 vmax = args.vmax,
                 rasterized=True)

    if args.invertcoastlines:
        m.drawcoastlines(linewidth = 0.25, color = 'white')
    else:
        m.drawcoastlines(linewidth = 0.25)

    if (args.colorbar):

        if (args.horizontal):
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes("bottom", size="5%", pad=0.05)
        
            if args.clabel is None:
                m.colorbar(p)
            else:
                m.colorbar(p, label = args.clabel)

        else:
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes("right", size="5%", pad=0.05)
        
            if args.clabel is None:
                m.colorbar(p)
            else:
                m.colorbar(p, label = args.clabel)

    if args.parallels:
        pminlat = (int(args.minlat/args.parallels) * args.parallels)
        if pminlat < args.minlat:
            pminlat = pminlat + args.parallels
            
        pmaxlat = (int(args.maxlat/args.parallels) * args.parallels)
        if pmaxlat > args.maxlat or pmaxlat == pminlat:
            pmaxlat = pmaxlat + args.parallels
        parallels = numpy.arange(pminlat, pmaxlat + args.parallels, args.parallels)
        #print 'Lat', pminlat, pmaxlat, parallels
        if args.invertcoastlines:
            m.drawparallels(parallels,labels=[1,0,0,1], color = 'w')
        else:
            m.drawparallels(parallels,labels=[1,0,0,1])

    if args.meridians:
        pminlon = (int(args.minlon/args.meridians) * args.meridians)
        if pminlon < args.minlon:
            pminlon = pminlon + args.meridians
            
        pmaxlon = (int(args.maxlon/args.meridians) * args.meridians)
        if pmaxlon < args.maxlon:
            pmaxlon = pmaxlon + args.meridians
        meridians = numpy.arange(pminlon, pmaxlon + args.meridians, args.meridians)
        #print 'Lon', pminlon, pmaxlon, meridians
        if args.invertcoastlines:
            m.drawmeridians(meridians,labels=[1,0,0,1], color = 'w')
        else:
            m.drawmeridians(meridians,labels=[1,0,0,1])

    if not args.points is None:

        tidepts, landpts, seapts = load_combined_points(args.points)

        print(args.point_filter)
        if args.point_filter == 0 or args.point_filter == 1:
            x = tidepts[:,0]
            y = tidepts[:,1]
            z = tidepts[:,2]

            m.scatter(x, y, c = z, latlon = True,
                      cmap = args.cmap,
                      vmin = args.vmin,
                      vmax = args.vmax,
                      edgecolor = 'black',
                      linewidth = 0.5)

        if args.point_filter == 0 or args.point_filter == 2:
            x = landpts[:,0]
            y = landpts[:,1]
            z = landpts[:,2]

            m.scatter(x, y, c = z, latlon = True,
                      cmap = args.cmap,
                      vmin = args.vmin,
                      vmax = args.vmax,
                      edgecolor = 'black',
                      linewidth = 0.5)
            
        if args.point_filter == 0 or args.point_filter == 3:
            x = seapts[:,0]
            y = seapts[:,1]
            z = seapts[:,2]

            m.scatter(x, y, c = z, latlon = True,
                      cmap = args.cmap,
                      vmin = args.vmin,
                      vmax = args.vmax,
                      edgecolor = 'black',
                      linewidth = 0.5)

    if args.pdf:
        fig.savefig(args.pdf, format = 'PDF')

    if args.png:
        fig.savefig(args.png, format = 'PNG', dpi = 300)

    if not args.noshow:
        P.show()
