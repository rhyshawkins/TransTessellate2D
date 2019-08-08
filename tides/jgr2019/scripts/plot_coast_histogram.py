
import argparse
import numpy
import matplotlib.pyplot as P

import math

def load_histogram(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    width, height, bins = map(int, lines[0].split())
    vmin, vmax = map(float, lines[1].split())

    counts = numpy.loadtxt(filename, skiprows = 2, dtype = 'int')
    #counts = numpy.array(map(lambda x: map(int, x.split()), lines[2:]), dtype = 'int')

    return width, height, bins, vmin, vmax, counts

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

def load_places(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    plons = []
    plats = []
    llons = []
    llats = []
    pnames = []
    
    for line in lines:
        t = line.split()

        plons.append(float(t[0]))
        plats.append(float(t[1]))
        llons.append(float(t[2]))
        llats.append(float(t[3]))

        pnames.append(' '.join(t[4:]).strip())

    return plons, plats, llons, llats, pnames

#
# This was copied from obspy code as the version on compute2 where this needed to
# be run has an outdated obspy.
#
WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563

def calc_vincenty_inverse(lat1, lon1, lat2, lon2, a=WGS84_A, f=WGS84_F):
    """
    Vincenty Inverse Solution of Geodesics on the Ellipsoid.

    Computes the distance between two geographic points on the WGS84
    ellipsoid and the forward and backward azimuths between these points.

    :param lat1: Latitude of point A in degrees (positive for northern,
        negative for southern hemisphere)
    :param lon1: Longitude of point A in degrees (positive for eastern,
        negative for western hemisphere)
    :param lat2: Latitude of point B in degrees (positive for northern,
        negative for southern hemisphere)
    :param lon2: Longitude of point B in degrees (positive for eastern,
        negative for western hemisphere)
    :param a: Radius of Earth in m. Uses the value for WGS84 by default.
    :param f: Flattening of Earth. Uses the value for WGS84 by default.
    :return: (Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees)
    :raises: This method may have no solution between two nearly antipodal
        points; an iteration limit traps this case and a ``StopIteration``
        exception will be raised.

    .. note::
        This code is based on an implementation incorporated in
        Matplotlib Basemap Toolkit 0.9.5 http://sourceforge.net/projects/\
matplotlib/files/matplotlib-toolkits/basemap-0.9.5/
        (matplotlib/toolkits/basemap/greatcircle.py)

        Algorithm from Geocentric Datum of Australia Technical Manual.

        * http://www.icsm.gov.au/gda/
        * http://www.icsm.gov.au/gda/gdatm/gdav2.3.pdf, pp. 15

        It states::

            Computations on the Ellipsoid

            There are a number of formulae that are available to calculate
            accurate geodetic positions, azimuths and distances on the
            ellipsoid.

            Vincenty's formulae (Vincenty, 1975) may be used for lines ranging
            from a few cm to nearly 20,000 km, with millimetre accuracy. The
            formulae have been extensively tested for the Australian region, by
            comparison with results from other formulae (Rainsford, 1955 &
            Sodano, 1965).

            * Inverse problem: azimuth and distance from known latitudes and
              longitudes
            * Direct problem: Latitude and longitude from known position,
              azimuth and distance.
    """
    # Check inputs
    if lat1 > 90 or lat1 < -90:
        msg = "Latitude of Point 1 out of bounds! (-90 <= lat1 <=90)"
        raise ValueError(msg)
    while lon1 > 180:
        lon1 -= 360
    while lon1 < -180:
        lon1 += 360
    if lat2 > 90 or lat2 < -90:
        msg = "Latitude of Point 2 out of bounds! (-90 <= lat2 <=90)"
        raise ValueError(msg)
    while lon2 > 180:
        lon2 -= 360
    while lon2 < -180:
        lon2 += 360

    b = a * (1 - f)  # semiminor axis

    if (abs(lat1 - lat2) < 1e-8) and (abs(lon1 - lon2) < 1e-8):
        return 0.0, 0.0, 0.0

    # convert latitudes and longitudes to radians:
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    TanU1 = (1 - f) * math.tan(lat1)
    TanU2 = (1 - f) * math.tan(lat2)

    U1 = math.atan(TanU1)
    U2 = math.atan(TanU2)

    dlon = lon2 - lon1
    last_dlon = -4000000.0  # an impossible value
    omega = dlon

    # Iterate until no significant change in dlon or iterlimit has been
    # reached (http://www.movable-type.co.uk/scripts/latlong-vincenty.html)
    iterlimit = 100
    try:
        while (last_dlon < -3000000.0 or dlon != 0 and
               abs((last_dlon - dlon) / dlon) > 1.0e-9):
            sqr_sin_sigma = pow(math.cos(U2) * math.sin(dlon), 2) + \
                pow((math.cos(U1) * math.sin(U2) - math.sin(U1) *
                     math.cos(U2) * math.cos(dlon)), 2)
            Sin_sigma = math.sqrt(sqr_sin_sigma)
            Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * \
                math.cos(U2) * math.cos(dlon)
            sigma = math.atan2(Sin_sigma, Cos_sigma)
            Sin_alpha = math.cos(U1) * math.cos(U2) * math.sin(dlon) / \
                math.sin(sigma)
            alpha = math.asin(Sin_alpha)
            Cos2sigma_m = math.cos(sigma) - \
                (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha), 2))
            C = (f / 16) * pow(math.cos(alpha), 2) * \
                (4 + f * (4 - 3 * pow(math.cos(alpha), 2)))
            last_dlon = dlon
            dlon = omega + (1 - C) * f * math.sin(alpha) * \
                (sigma + C * math.sin(sigma) *
                    (Cos2sigma_m + C * math.cos(sigma) *
                        (-1 + 2 * pow(Cos2sigma_m, 2))))

            u2 = pow(math.cos(alpha), 2) * (a * a - b * b) / (b * b)
            A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
            B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
            delta_sigma = B * Sin_sigma * \
                (Cos2sigma_m + (B / 4) *
                    (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) - (B / 6) *
                        Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) *
                        (-3 + 4 * pow(Cos2sigma_m, 2))))

            dist = b * A * (sigma - delta_sigma)
            alpha12 = math.atan2(
                (math.cos(U2) * math.sin(dlon)),
                (math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) *
                 math.cos(dlon)))
            alpha21 = math.atan2(
                (math.cos(U1) * math.sin(dlon)),
                (-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) *
                 math.cos(dlon)))
            iterlimit -= 1
            if iterlimit < 0:
                # iteration limit reached
                print('Iteration', lon1, lat1, lon2, lat2)
                raise StopIteration
    except ValueError:
        # usually "math domain error"
        print('ValueError')
        raise StopIteration

    if alpha12 < 0.0:
        alpha12 = alpha12 + (2.0 * math.pi)
    if alpha12 > (2.0 * math.pi):
        alpha12 = alpha12 - (2.0 * math.pi)

    alpha21 = alpha21 + math.pi

    if alpha21 < 0.0:
        alpha21 = alpha21 + (2.0 * math.pi)
    if alpha21 > (2.0 * math.pi):
        alpha21 = alpha21 - (2.0 * math.pi)

    # convert to degrees:
    alpha12 = alpha12 * 360 / (2.0 * math.pi)
    alpha21 = alpha21 * 360 / (2.0 * math.pi)

    return dist, alpha12, alpha21


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--coast-points', type = str, required = True, help = 'Coastal Points')
    parser.add_argument('-C', '--coast-places', type = str, default = None, help = 'Coastal places')
                        
    parser.add_argument('-H', '--histogram', type = str, required = True, help = 'Histogram')

    parser.add_argument('-d', '--data', type = str, default = None, help = 'Data')

    parser.add_argument('--minlon', type = float, default = -10.0, help = 'Minlon')
    parser.add_argument('--maxlon', type = float, default = 10.0, help = 'Maxlon')
    parser.add_argument('--minlat', type = float, default = -10.0, help = 'Minlat')
    parser.add_argument('--maxlat', type = float, default = 10.0, help = 'Maxlat')

    parser.add_argument('--true-sea', type = str, default = None, help = 'True sea image')
    parser.add_argument('--true-land', type = str, default = None, help = 'True land image')

    parser.add_argument('--mean-tide', type = str, default = None, help = 'Mean tide image')
    
    parser.add_argument('--roll', type = int, default = 0, help = 'Roll points')                    

    parser.add_argument('--log', action = 'store_true', default = False, help = 'Log histogram')
    parser.add_argument('--log-threshold', type = float, default = 1.0e-6, help = 'Min probability')
    
    parser.add_argument('--pdf', type = str, default = None, help = 'PDF output')
    parser.add_argument('--width', type = float, default = 8.0, help = 'Plot width')
    parser.add_argument('--height', type = float, default = 6.0, help = 'Plot height')

    parser.add_argument('--colorbar', action = 'store_true', default = False, help = 'Show colour bar')

    parser.add_argument('--vmin', type = float, default = None, help = 'V min')
    parser.add_argument('--vmax', type = float, default = None, help = 'V max')

    parser.add_argument('--cmap', type = str, default = 'jet', help = 'Color map')

    parser.add_argument('--absolute', action = 'store_true', default = False, help = 'Abs. sea level')
    
    args = parser.parse_args()

    points = numpy.loadtxt(args.coast_points)

    npoints, _ = points.shape

    if (args.roll != 0):

        # We keep the distance column the same (approximately evenly spaced points
        # so it should be reasonably accurate.
        points[:,0] = numpy.roll(points[:,0], args.roll)
        points[:,1] = numpy.roll(points[:,1], args.roll)


    width, height, bins, vmin, vmax, counts = load_histogram(args.histogram)

    image = numpy.zeros((bins, npoints))

    x = numpy.linspace(0.0, points[-1,2], npoints) # Approximation

    for i in range(npoints):

        xi = int((points[i, 0] - args.minlon)/(args.maxlon - args.minlon) * float(width))
        if (xi < 0):
            xi = 0
        if (xi >= width):
            xi = width - 1

        yi = int((points[i, 1] - args.minlat)/(args.maxlat - args.minlat) * float(height))
        if (yi < 0):
            yi = 0
        if (yi >= height):
            yi = height - 1

        t = numpy.array(counts[yi * width + xi,:], dtype = 'float', copy = True)
        t = (t/numpy.sum(t)) / (vmax - vmin)

        if args.log:
            valid_indices = numpy.where(t >= args.log_threshold)[0]
            inv_indices = numpy.where(t < args.log_threshold)[0]
            t[valid_indices] = numpy.log10(t[valid_indices])
            t[inv_indices] = numpy.log10(args.log_threshold)
            image[:,i] = t
        else:
            image[:,i] = t

    fig, ax = P.subplots()
    fig.set_tight_layout(True)
    fig.set_size_inches(args.width, args.height)

    imin = numpy.min(image)
    imax = numpy.max(image)
    if not (args.vmin is None):
        imin = args.vmin
    if not (args.vmax is None):
        imax = args.vmax
    
    g = ax.imshow(image, vmin = imin, vmax = imax,
                  extent = [x[0], x[-1], vmin, vmax],
                  aspect = 'auto',
                  origin = 'lower',
                  cmap = args.cmap)

    #
    # Plotting data
    #
    if not (args.data is None):

        tides, _, _ = load_combined_points(args.data)

        ndata, _ = tides.shape

        dx = numpy.zeros((ndata,))

        for i in range(ndata):

            mindist = 1e30
            mindistj = -1
            # Find nearest point

            for j in range(npoints):

                dist, _, _ = calc_vincenty_inverse(tides[i, 1], tides[i, 0], points[j, 1], points[j, 0])

                dist = dist/1.0e3
                
                if dist < mindist:
                    mindist = dist
                    mindistj = j

            dx[i] = points[mindistj, 2]

#            if (dist > 2000.0):
#                print dist
#                dx[i] = -1
#            else:


        indices = numpy.where(dx >= 0)[0]

        ax.errorbar(dx[indices], tides[indices,2], yerr = tides[indices,3], linestyle = 'none', capsize = 5)

    #
    # Plotting true curve
    #
    if not ((args.true_sea is None) or (args.true_land is None)):

        tsea = numpy.loadtxt(args.true_sea)
        tland = numpy.loadtxt(args.true_land)

        if tsea.shape != tland.shape:
            raise Exception('True image size mismatch')
        
        theight, twidth = tsea.shape
        
        ttide = numpy.zeros((npoints, ))
        
        for i in range(npoints):

            xi = int((points[i, 0] - args.minlon)/(args.maxlon - args.minlon) * float(twidth))
            if (xi < 0):
                xi = 0
            if (xi >= width):
                xi = width - 1
                    
            yi = int((points[i, 1] - args.minlat)/(args.maxlat - args.minlat) * float(theight))
            if (yi < 0):
                yi = 0
            if (yi >= height):
                yi = height - 1
                    

            ttide[i] = tsea[yi, xi] - tland[yi, xi]

        ax.plot(points[:,2], ttide, 'w:')

    #
    # Plotting mean curve
    #
    if not (args.mean_tide is None):

        t = numpy.loadtxt(args.mean_tide)

        theight, twidth = t.shape
        
        ttide = numpy.zeros((npoints, ))
        
        for i in range(npoints):

            xi = int((points[i, 0] - args.minlon)/(args.maxlon - args.minlon) * float(twidth))
            if (xi < 0):
                xi = 0
            if (xi >= width):
                xi = width - 1
                    
            yi = int((points[i, 1] - args.minlat)/(args.maxlat - args.minlat) * float(theight))
            if (yi < 0):
                yi = 0
            if (yi >= height):
                yi = height - 1
                    

            ttide[i] = t[yi, xi]

        ax.plot(points[:,2], ttide, 'g--')

    #
    # Place names
    #
    if not (args.coast_places is None):


        plons, plats, llons, llats, placenames = load_places(args.coast_places)
        
            
        N = len(placenames)
        for i in range(N):
            lon = plons[i]
            lat = plats[i]

            mindist = 1e30
            mindistj = -1
            # Find nearest point

            for j in range(npoints):
                dlon = points[j, 0] - lon
                dlat = points[j, 1] - lat
                dist = numpy.sqrt(dlon*dlon + dlat*dlat)

                if dist < mindist:
                    mindist = dist
                    mindistj = j

            x = points[mindistj, 2]

            placecolor = 'yellow'
            ax.axvline(x = x, linestyle = 'solid', color = placecolor, linewidth = 0.5)

            if i == (N - 1):
                ax.text(x, vmin + 0.01 * (vmax - vmin),
                        placenames[i], va = 'bottom', ha = 'right', rotation = 90, color = placecolor, size = 8)
            else:
                ax.text(x, vmin + 0.01 * (vmax - vmin),
                        placenames[i], va = 'bottom', ha = 'left', rotation = 90, color = placecolor, size = 8)
            
    ax.set_xlim(0, numpy.max(points[:,2]))

    ax.set_ylim(vmin, vmax)

    if args.absolute:
        ax.set_ylabel('Abs. Sea Level Rate (mm/year)')
    else:
        ax.set_ylabel('Rel. Sea Level Rate (mm/year)')
        
    ax.set_xlabel('Distance along coast (km)')
    
    if args.colorbar:
        P.colorbar(g, ax = ax)

    if args.pdf is None:
        P.show()
    else:
        fig.savefig(args.pdf, format = 'PDF')
        
        
        
