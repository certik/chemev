"""Bezier curve calculation and interpolation.
2006 Ondrej Certik & Henry Ferguson, STScI

Usage:

    import bezier
    x,y=bezier.curve((
        bezier.point((0,1),(1,0)),
        bezier.point((1,0),(-0.5,0.5)),
        bezier.point((1.5,0),(-0.5,0.5),(0,0.5)),
        bezier.point((2,1),(-0.5,0))
        ))
    import pylab
    pylab.plot(x,y)
    pylab.show()

bezier.point() class' first parameter is the coordinate of the point,
second is the in tangent (as a vector, relative to the point) and 
third is the out tangent. if the outtangent is omitted, it is taken as
"-intangent". for the first and last point, you must just specify the
intangent.

you can specify the level of the bezier curve by the optional parameter to the
curve function:

def curve(points,level=6):


"""

import numarray

def midpoint((x1, y1), (x2, y2)):
    return ((x1+x2)/2., (y1+y2)/2.)

class bezier:
    def __init__(self,a,b,max_level=6):
        self.x0 = a[0]
        self.y0 = a[1] 
        self.x1 = b[0]
        self.y1 = b[1]
        self.coordsp1 = []
        self.coordsp4 = []
        self.max_level = max_level

    def draw_curve(self,P1, P2, P3, P4, level=1):
        if level == self.max_level:
            self.coordsp1 += [P1]
            self.coordsp4 += [P4]
        else:
            L1 = P1
            L2 = midpoint(P1, P2)
            H  = midpoint(P2, P3)
            R3 = midpoint(P3, P4)
            R4 = P4
            L3 = midpoint(L2, H)
            R2 = midpoint(R3, H)
            L4 = midpoint(L3, R2)
            R1 = L4
            self.draw_curve(L1, L2, L3, L4, level+1)
            self.draw_curve(R1, R2, R3, R4, level+1)

    def bezier(self,a,b):
        self.coordsp1 = []
        self.coordsp4 = []
        self.draw_curve((self.x0,self.y0),(a[0],a[1]),(b[0],b[1]),(self.x1,self.y1))
        xy = numarray.array(self.coordsp1+self.coordsp4[-1:])
        self.x = xy[:,0]
        self.y = xy[:,1]
        return self.x,self.y

def curve(points,level=6):
    points[0].middle=False
    points[-1].middle=False
    pold=points[0]
    x=[pold.pos[0]]; y=[pold.pos[1]]
    for pnew in points[1:]:
        x0,y0=bezier(pold.pos,pnew.pos,level).bezier(
                pold.getmiddle2(),pnew.getmiddle1())
        x.extend(x0[1:])
        y.extend(y0[1:])
        pold=pnew
    return numarray.array(x),numarray.array(y)

class point:
    def __init__(self,pos,intangent=None,outtangent=None):
        self.pos=pos
        self.intangent=intangent
        self.outtangent=outtangent
        self.middle=True

    def add(self,vec):
        return (self.pos[0]+vec[0],self.pos[1]+vec[1])

    def getfirstlast(self):
        if self.intangent!=None:
            assert self.outtangent==None
            return self.add(self.intangent)
        elif self.outtangent!=None:
            assert self.intangent==None
            return self.add(self.outtangent)
        else: 
            raise "need to specify exactly one tangent"

    def getmiddle1(self):
        if not self.middle: 
            return self.getfirstlast()
        if self.intangent!=None:
            return self.add(self.intangent)
        elif self.outtangent!=None:
            return self.add((-self.outtangent[0],-self.outtangent[1]))
        else:
            raise "need to specify exactly one tangent"

    def getmiddle2(self):
        if not self.middle: 
            return self.getfirstlast()
        if self.outtangent!=None:
            return self.add(self.outtangent)
        elif self.intangent!=None:
            return self.add((-self.intangent[0],-self.intangent[1]))
        else:
            raise "need to specify exactly one tangent"



def interpolate(points,t):
    """points=(x,y), list of x's and y's. "t" ... list of x, for which
    we want to evaluate the dependence "y". returns t,yy, where yy is the
    "y" evaluated at the grid "t".
    """
    xx,yy=points
    x=t

    t0 = numarray.compress(x<=min(xx), x)
    t1 = numarray.compress( (x>min(xx)) & (x<max(xx)), x )
    t2 = numarray.compress(x>=max(xx), x)
    slope0 = (yy[1]-yy[0])/(xx[1]-xx[0])
    slope2 = (yy[-1]-yy[-2])/(xx[-1]-xx[-2])
    indices = numarray.searchsorted(xx,t1)
    x0 = xx[indices-1]
    x1 = xx[indices]
    y0 = yy[indices-1]
    y1 = yy[indices]
    slope = (y1-y0)/(x1-x0)

    y1 = slope*(t1-x0)+y0
    y0 = slope0*(t0-xx[0])+yy[0]    
    y2 = slope2*(t2-xx[-1])+yy[-1]    

    y = numarray.concatenate((y0,y1,y2))
    return t,y
