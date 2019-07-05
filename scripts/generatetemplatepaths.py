
import math
import sys
import optparse
import random

import tomosynthetic

EARTH_RADIUS = 6371.0

def vrot(v, theta):
    #c -s
    #s  c

    x, y = v
    rad = theta * math.pi/180.0
    
    xp = math.cos(rad)*x - math.sin(rad)*y
    yp = math.sin(rad)*x + math.cos(rad)*y

    return (xp, yp)

def bezier(p0, p1, p2, t):
    omt = 1.0 - t
    return omt*omt*p0 + 2.0*omt*t*p1 + t*t*p2

def generate_path(lon1, lat1, r1, lon2, lat2, r2, delta, v):

    clon = (lon1 + lon2)/2.0
    clat = (lat1 + lat2)/2.0

    dlon = (lon2 - lon1)
    dlat = (lat2 - lat1)

    n = int(math.ceil(math.sqrt(dlon*dlon + dlat*dlat)/delta))
    # Ensure we have at least a few points
    if n < 5: 
        n = 5
        
    if (v == 0.0):
        #
        # Straight line
        #

        p = map(lambda i: (lon1 + dlon*float(i)/float(n - 1),
                           lat1 + dlat*float(i)/float(n - 1),
                           r1 + (r2 - r1)*float(i)/float(n - 1)), range(n))

    else:

        if (v < 0.0):
            v = -v
            nlon, nlat = vrot((dlon, dlat), -90.0)
        else:
            nlon, nlat = vrot((dlon, dlat), 90.0)

        nlen = math.sqrt(nlon*nlon + nlat*nlat)

        plon = clon + v*nlon/nlen
        plat = clat + v*nlat/nlen

        p = map(lambda i: (bezier(lon1, plon, lon2, float(i)/float(n - 1)),
                           bezier(lat1, plat, lat2, float(i)/float(n - 1)),
                           r1 + (r2 - r1)*float(i)/float(n - 1)), range(n))
        
    return p
        


def sphtocart(lon, lat, r):
    lon = lon * math.pi/180.0
    lat = lat * math.pi/180.0

    x = r * math.cos(lon) * math.cos(lat)
    y = r * math.sin(lon) * math.cos(lat)
    z = r * math.sin(lat)

    return (x, y, z)
    
def veclen(x, y, z):
    return math.sqrt(x*x + y*y + z*z)

def cartesian_dist(lon1, lat1, r1, lon2, lat2, r2):

    x1, y1, z1 = sphtocart(lon1, lat1, r1)
    x2, y2, z2 = sphtocart(lon2, lat2, r2)

    return veclen(x2 - x1, y2 - y1, z2 - z1)
    
def trace_path(path, velocity_fn):

    dist = [0.0] * len(path)
    
    # First pass, compute distances
    for i in range(1, len(path)):
        
        lon1, lat1, r1 = path[i - 1]
        lon2, lat2, r2 = path[i]

        d = cartesian_dist(lon1, lat1, r1, lon2, lat2, r2)

        dist[i - 1] = dist[i - 1] + d/2.0
        dist[i] = d/2.0
    
    # 2nd pass, integrate travel time
    t = 0.0
    total_dist = 0.0
    mean_v = 0.0
    for i in range(len(path)):

        lon, lat, r = path[i]
        v = velocity_fn(lon, lat, r)
        
        total_dist = total_dist + dist[i]
        mean_v = mean_v + v
        t = t + dist[i]/v

    mean_v = mean_v/float(len(path))
    print '      ', t, total_dist, total_dist/t, mean_v, total_dist/mean_v
    return t

if __name__ == '__main__':

    parser = optparse.OptionParser()

    parser.add_option('--minlon', dest='minlon', type='float', default = -10.0, help='Lower longitude bound')
    parser.add_option('--maxlon', dest='maxlon', type='float', default = 10.0, help='Upper longitude bound')
    parser.add_option('--minlat', dest='minlat', type='float', default = -10.0, help='Lower latitude bound')
    parser.add_option('--maxlat', dest='maxlat', type='float', default = 10.0, help='Upper latitude bound')

    parser.add_option('--radius', dest='R', type='float', default = EARTH_RADIUS, help='Min radius')

    parser.add_option('-v', '--variance', dest='variance', type='float', default = 1.0, help='Variance of curves (0 = straight lines.')
    parser.add_option('-d', '--delta', dest='delta', type='float', default = 0.2, help='Path step size')

    parser.add_option('-p', '--paths', dest='paths', type='int', default = 100, help='No. of observations to create')
    parser.add_option('-n', '--noise', dest='noise', type='float', default = 1.0, help='Std. deviation of noise in observed time.')

    parser.add_option('-o', '--observations', dest='observed', default = None, help='Observations file to write')
    parser.add_option('--real', dest='real', default = None, help='Exact Observations file to write')

    parser.add_option('-t', '--tomography', dest='tomography', default = tomosynthetic.DEFAULT_MODEL, help='Tomography model to use')

    parser.add_option('-l', '--list-tomography', dest='listtomography', action='store_true', help='List available tomography models')

    parser.add_option('--data', dest='data', default = None, help = 'Write data to file.')

    parser.add_option('-N', '--enable-noise', dest='enablenoise', default = False, action='store_true', help='Add noise to observations.')

    options, args = parser.parse_args()

    if (options.listtomography):
        print 'Default model:\n    "%s"' % tomosynthetic.DEFAULT_MODEL
        print 'Available models:'
        for k in tomosynthetic.fmap.keys():
            if len(k):
                print '    "%s"' % k
        sys.exit(0)

    if (options.observed == None):
        print 'Missing observed file option.'
        sys.exit(-1)

    if (not tomosynthetic.fmap.has_key(options.tomography)):
        print 'Invalid tomography model.'
        sys.exit(-1)

    #
    # Create the random paths and compute traveltimes
    #
    sources = []
    receivers = []
    obs = []
    f = open(options.observed, 'w')
    f.write('%d\n' % options.paths)
    for i in range(options.paths):
        slon = random.uniform(options.minlon, options.maxlon)
        slat = random.uniform(options.minlat, options.maxlat)
        sources.append((slon, slat))
        rlon = random.uniform(options.minlon, options.maxlon)
        rlat = random.uniform(options.minlat, options.maxlat)
        receivers.append((rlon, rlat))

        p = generate_path(slon, slat, options.R, rlon, rlat, options.R, options.delta, options.variance)

        N = 0.0
        if options.enablenoise:
            N = random.normalvariate(0.0, options.noise)

        t = trace_path(p, tomosynthetic.fmap[options.tomography])
        
        obs.append((p, t))

        f.write('%d %f %f\n' % (len(p), t + N, options.noise))
        for (lon, lat, r) in p:
            f.write('%.10f %.10f\n' % (lon, lat))
    
    f.close()

    if options.data:
        f = open(options.data, 'w')

        xsamples = 100
        ysamples = 100
        for ilat in range(ysamples):
            lat = float(ilat)/float(ysamples - 1) * (options.maxlat - options.minlat) + options.minlat
            for ilon in range(xsamples):
                lon = float(ilon)/float(xsamples - 1) * (options.maxlon - options.minlon) + options.minlon

                f.write('%f ' % tomosynthetic.fmap[options.tomography](lon, lat, tomodata.EARTH_RADIUS))

            f.write('\n')

        f.close()
        
            
    


        

    
