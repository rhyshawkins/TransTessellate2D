
import sys

import argparse
import numpy

def sdf_filter(sdf, xmin, xmax, ymin, ymax, boundary, f):

    points = []

    r, c = sdf.shape

    for j in range(r):
        for i in range(c):

            if (f(sdf[j, i])):

                lon = xmin + (float(i)+0.5)/float(c) * (xmax - xmin)
                lat = ymin + (float(j)+0.5)/float(r) * (ymax - ymin)

                if (lon > (xmin + boundary) and lon < (xmax - boundary) and
                    lat > (ymin + boundary) and lat < (ymax - boundary)):
                    
                    points.append((lon, lat))

    return points
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--sdf', type = str, required = True, help = 'SDF file')
    
    parser.add_argument('-L', '--land', type = int, default = 50, help = 'No. land observations')

    parser.add_argument('-S', '--sea', type = int, default = 50, help = 'No. sea observations')

    parser.add_argument('-T', '--tide', type = int, default = 50, help = 'No. tide observations')
    
    parser.add_argument('--xmin', type = float, default = -10.0, help = 'Min x')
    parser.add_argument('--ymin', type = float, default = -10.0, help = 'Min y')
    parser.add_argument('--xmax', type = float, default = 10.0, help = 'Max x')
    parser.add_argument('--ymax', type = float, default = 10.0, help = 'Max y')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('--threshold-sea', type = float, default = 10.0, help = 'Distance from shore for altimetry')
    parser.add_argument('--threshold-land', type = float, default = 5.0, help = 'Distance from shore for tide gauges')

    parser.add_argument('--tide-boundary', type = float, default = 0.0, help = 'Boundary buffer for tide gauges (degrees)')
    
    parser.add_argument('--land-boundary', type = float, default = 0.0, help = 'Boundary buffer for GPS (degrees)')

    args = parser.parse_args()

    sdf = numpy.loadtxt(args.sdf)

    possible_sea = sdf_filter(sdf,
                              args.xmin,
                              args.xmax,
                              args.ymin,
                              args.ymax,
                              0.0,
                              lambda x: x > args.threshold_sea)

    possible_tide = sdf_filter(sdf,
                               args.xmin,
                               args.xmax,
                               args.ymin,
                               args.ymax,
                               args.tide_boundary,
                               lambda x: x < 0.0 and x > -args.threshold_land)

    possible_land = sdf_filter(sdf,
                               args.xmin,
                               args.xmax,
                               args.ymin,
                               args.ymax,
                               args.land_boundary,
                               lambda x: x < -args.threshold_land)

    if (args.land > len(possible_land)):
        print 'error: only %d land points available and %d requested' % (len(possible_land), args.land)
        sys.exit(-1)

    if (args.sea > len(possible_sea)):
        print 'error: only %d sea points available and %d requested' % (len(possible_sea), args.sea)
        sys.exit(-1)

    if (args.tide > len(possible_tide)):
        print 'error: only %d tide points available and %d requested' % (len(possible_tide), args.tide)
        sys.exit(-1)

    numpy.random.seed(983)
    
    numpy.random.shuffle(possible_sea)
    numpy.random.shuffle(possible_land)
    numpy.random.shuffle(possible_tide)
        
    f = open(args.output, 'w')

    nobs = args.land + args.sea + args.tide
    
    f.write('%d\n' % nobs)

    for lon, lat in possible_tide[:args.tide]:

        f.write('%15.9f %15.9f 1 0.0 1.0\n' % (lon, lat))
        
    for lon, lat in possible_land[:args.land]:

        f.write('%15.9f %15.9f 2 0.0 1.0\n' % (lon, lat))

    for lon, lat in possible_sea[:args.sea]:

        f.write('%15.9f %15.9f 3 0.0 1.0\n' % (lon, lat))
        
    f.close()
    

    
    
