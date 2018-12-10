
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-L', '--land', type = int, default = 50, help = 'No. land observations')

    parser.add_argument('-S', '--sea', type = int, default = 50, help = 'No. sea observations')

    parser.add_argument('-T', '--tide', type = int, default = 50, help = 'No. tide observations')
    
    parser.add_argument('--xmin', type = float, default = -10.0, help = 'Min x')
    parser.add_argument('--ymin', type = float, default = -10.0, help = 'Min y')
    parser.add_argument('--xmax', type = float, default = 10.0, help = 'Max x')
    parser.add_argument('--ymax', type = float, default = 10.0, help = 'Max y')

    parser.add_argument('--radius', type = float, default = 5.0, help = 'Island radius')

    parser.add_argument('--half-tide', action = 'store_true', default = False, help = 'Single half of island has tide gauges')
    parser.add_argument('--half-land', action = 'store_true', default = False, help = 'Single half of island has land gps')
    
    
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    args = parser.parse_args()

    f = open(args.output, 'w')

    nobs = args.land + args.sea + args.tide
    
    f.write('%d\n' % nobs)

    cx = (args.xmax + args.xmin)/2.0
    cy = (args.ymax + args.ymin)/2.0

    thetarange = numpy.pi * 2.0
    if args.half_tide:
        thetarange = numpy.pi
        
    for i in range(args.tide):

        r = args.radius
        theta = numpy.random.uniform() * thetarange

        x = -r * numpy.sin(theta)
        y = r * numpy.cos(theta)
        
        f.write('%15.9f %15.9f 1 0.0 1.0\n' % (cx + x, cy + y))
        
    thetarange = numpy.pi * 2.0
    if args.half_land:
        thetarange = numpy.pi

    for i in range(args.land):

        r = numpy.random.uniform() * args.radius
        theta = numpy.random.uniform() * thetarange
        

        x = -r * numpy.sin(theta)
        y = r * numpy.cos(theta)
        
        f.write('%15.9f %15.9f 2 0.0 1.0\n' % (cx + x, cy + y))

    for i in range(args.sea):

        r = 0.0
        while r < args.radius:

            x = args.xmin + (args.xmax - args.xmin) * numpy.random.uniform()
            y = args.ymin + (args.ymax - args.ymin) * numpy.random.uniform()

            dx = x - cx
            dy = y - cy
            
            r = numpy.sqrt(dx*dx + dy*dy)

        f.write('%15.9f %15.9f 3 0.0 1.0\n' % (x, y))
        
    f.close()
    

    
    
