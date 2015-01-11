import body, hipparcos, constellation
import math
import xml.dom.minidom as dom
import itertools
import datetime

def conv(i):
    if isinstance(i, float):
        return '%.1f'%i
    return str(i)

class Chart:
    # These are global class declarations.
    inches = 12
    radius = 1000*inches/12
    daytick_space = 25*inches/12
    daytick_utchour = 4
    daytick_monthsize = 48*inches/12
    padding = 10
    doc = dom.Document()
    hip_stars = None
    # What is the size of a zero magnitude star?
    zerosize = math.sqrt(inches*6)

    def __init__(self, loc, begin, end):
        self.init_svg()
        self.init_orientation(loc)
        self.read_stars()
        # Produce the dayticks.
        self.init_dayticks(loc, begin, end)
        # Produce the gridlines.
        self.gridlines(15)

    def init_dayticks(self, loc, begin, end):
        # Produce the outermost circles.
        c = self.doc.createElement('circle')
        attr = {'cx':0, 'cy':0, 'r':self.radius, 'stroke-width':4}
        for k,v in attr.items(): c.setAttribute(k,conv(v))
        self.centered.appendChild(c)
        tickgroup = self.make_element(self.centered, 'g', ('stroke-width',2))
        self.make_element(tickgroup, 'circle', (
            'r',self.radius), ('stroke-width',4))
        self.make_element(tickgroup, 'circle', (
            'r',self.radius+self.daytick_space), ('stroke-width',4))
        # Get the appropriate year.
        year = datetime.datetime.utcnow().year
        # Now, for each day of this year...
        dt = datetime.datetime(year, 1, 1, self.daytick_utchour)
        timedelta = datetime.timedelta(1)
        while dt.year==year:
            dn = body.Convert.date2day_number(dt)
            sidereal = loc.sidereal(dn)
            x,y = self.ra2xy(sidereal)
            line = self.make_element(tickgroup, 'line', (
                'x1',x*self.radius), ('y1',y*self.radius), (
                'x2',x*(self.radius+self.daytick_space)), (
                'y2',y*(self.radius+self.daytick_space)))
            if dt.day==1:
                line.setAttribute('stroke-width', '6')
            dt += timedelta
        # For each month, create the label.
        monthgroup = self.make_element(
            self.centered, 'g', ('stroke', 'none'),
            ('font-size', '%dpt'%(self.daytick_monthsize-20)),
            ('fill', 'black'), ('text-anchor', 'middle'))
        for m in xrange(1,13):
            dt = datetime.datetime(year, m, 1, self.daytick_utchour)
            sr = loc.sidereal(body.Convert.date2day_number(dt))
            dtnext = datetime.datetime(year, (m%12)+1, 1,self.daytick_utchour)
            sr_next = loc.sidereal(body.Convert.date2day_number(dtnext))
            (x1,y1), (x2,y2) = self.ra2xy(sr), self.ra2xy(sr_next)
            radius = self.radius+self.daytick_space+4
            x1,y1,x2,y2 = (i*radius for i in (x1,y1,x2,y2))
            self.make_element(
                self.defs, 'path',
                ('d', 'M%d,%d A%d,%d 0 0,1 %d,%d'%(
                x1,y1,radius+10,radius+10,x2,y2)), (
                'id', 'montharc-%d'%m))

            #'d', 'M%d,%d A%d,%d 0 1,1 %d,%d' % (
            #x2,y2,radius,radius,x1,y1)), (
            t = self.make_element(monthgroup, 'text')
            tp = self.make_element(t, 'textPath', (
                'xlink:href','#montharc-%d'%m), ('startOffset', '50%'))
            tp.appendChild(self.doc.createTextNode(dt.strftime('%B').upper()))

    @classmethod
    def make_element(self, parent, name, *attrs):
        """Create an element with a given parent, tag name, and attributes.

        This will create a given XML element with a tag name given as
        a string.  Attributes are set as optional arguments after the
        parent and name, where each argument is represented as a
        two-element tuple of name (a string) and value.  If the parent
        argument is set to None (or otherwise evalutes to false) then
        the created node shall not be added to any parent."""
        e = self.doc.createElement(name)
        for k,v in attrs: e.setAttribute(k,conv(v))
        if parent: parent.appendChild(e)
        return e

    @classmethod
    def make_textelement(self, parent, text, *attrs):
        t = self.make_element(parent, 'text', *attrs)
        t.appendChild(self.doc.createTextNode(text))
        return t

    def init_svg(self):
        """An internal method for initializing the SVG.

        At the end of this call, the attributes 'svg' (the SVG dom
        element) and 'centered' (the G element with the center of the
        chart at coordinates (0,0)) will be defined.  The link to the
        defs where referenced paths are declared shall also be created
        as the 'defs' member.  The image shall be empty in terms of
        any visual content."""
        self.svg = self.doc.createElement('svg')
        halfwidth = self.radius+self.daytick_space+self.daytick_monthsize+\
                    self.padding
        dimension = 2*halfwidth
        attr = {'xmlns':'http://www.w3.org/2000/svg', 'version':'1.1',
                'xmlns:xlink':'http://www.w3.org/1999/xlink',
                'viewBox':'0 0 %d %d'%(dimension,dimension),
                'height':'%din'%self.inches, 'width':'%din'%self.inches, 
                'preserveAspectRatio':'xMinYMid meet',
                'stroke':'black', 'fill':'none',
                'font-family':'Arial', 'font-size':10}
        for k,v in attr.items(): self.svg.setAttribute(k,conv(v))
        # Create the clipping path for the interior region of the chart.
        self.defs = self.make_element(self.svg, 'defs')
        clip = self.make_element(
            self.defs, 'clipPath', ('id', 'innerClipPath'))
        self.make_element(
            clip, 'circle', ('cx',0), ('cy',0), ('r',self.radius))
        # Make 0,0 the center of the circle.
        self.centered = self.doc.createElement('g')
        self.centered.setAttribute('transform','translate(%d,%d)'%(
            2*(halfwidth,)))
        self.svg.appendChild(self.centered)

    def init_orientation(self, loc):
        """Set the declination in radians of regions of the chart.
        
        At the end of this call, these attributes shall be set:

        northern - A boolean value indicating whether the chart is for
        the northern hemisphere.
        
        inner_dec - The declination of the innermost portion of the
        chart.  This corresponds to the celestial north or south
        (depending on hemisphere) and is correspondingly only pi/2 or
        -pi/2 (depending on pole).

        outer_dec - The declination of the outermost portion of the
        chart.  This corresponds to the minimum (or maximum)
        declination visible in the sky from the indicated location,
        which can be seen by facing directly south (or north) if we
        are in the north (or south) hemisphere.

        Note that these latter two quantities are in radians."""
        # The chart is arranged in a circle, with the innermost circle
        # that of the visible celestial pole.
        self.northern = loc.la >= 0
        if self.northern:
            # If we are in the northern hemisphere, the innermost portion
            # of the chart circle should be the celestial north pole, and
            # the outermost should correspond to the declination at the
            # south horizon.
            self.inner_dec, self.outer_dec = math.pi/2.0, loc.la-math.pi/2.0
        else:
            # If we are in the southern hemisphere, the innermost portion
            # of the chart circle should be the celestial south pole, and
            # the outermost should correspond to the declination at the
            # northern horizon.
            self.inner_dec, self.outer_dec = -math.pi/2.0, loc.la+math.pi/2.0

    @classmethod
    def read_stars(self):
        """Reads in all the stars.

        After this class method is called, the stars should be in the
        hip_stars attribute of the class.  If this attribute is
        already defined with a non-trivial argument, this method shall
        do nothing."""
        if self.hip_stars: return
        all_stars = list(hipparcos.stars())
        self.hip_stars = [None]*(max(s[0] for s in all_stars)+1)
        for s in all_stars: self.hip_stars[s[0]] = s

    def dec2radius(self, dec):
        """For a declination, get the corresponding distance from center.

        For a declination given in radians, this will produce the ."""
        return self.radius*(self.inner_dec-dec)/(self.inner_dec-self.outer_dec)

    def dec2radius(self, dec):
        if self.northern:
            inner, outer = math.pi/2-dec, math.pi/2-self.outer_dec
        else:
            inner, outer = math.pi/2+dec, math.pi/2+self.outer_dec
        return self.radius*math.sqrt((
            1.0-math.cos(inner))/(1.0-math.cos(outer)))

    def ra2xy(self, ra):
        """Given a RA, give an x, y pair ratio of offset from center.

        For example, returning .6, -.8 would indicate that the region
        of the graphic in increments right 3 and up 4 will all
        correspond to the RA given as an argument."""
        return -math.sin(ra), math.cos(ra)

    def radec2xy(self, ra, dec):
        r = self.dec2radius(dec)
        return -r*math.sin(ra), r*math.cos(ra)

    def gridlines(self, spacing, stroke=.5):
        """Produce the concentric circles."""
        # Produce the g element for the gridlines.
        gridline_g = self.doc.createElement('g')
        gridline_g.setAttribute('stroke-width', '%g'%stroke)
        self.centered.appendChild(gridline_g)
        # Produce the g element for the labels.
        labels_g = self.doc.createElement('g')
        attr = {'fill':'black', 'stroke':'none'}
        for k,v in attr.items(): labels_g.setAttribute(k,conv(v))
        self.centered.appendChild(labels_g)
        # Now draw them.
        innermost = 90-spacing if self.northern else -90+spacing
        steps = xrange(innermost, int(math.degrees(self.outer_dec)), -spacing)
        for circle in steps:
            radius = self.dec2radius(math.radians(circle))
            # First create the circle.
            c = self.doc.createElement('circle')
            attr = {'cx':0, 'cy':0, 'r':radius}
            for k,v in attr.items(): c.setAttribute(k,conv(v))
            gridline_g.appendChild(c)
            # Then create the label.
            t = self.make_textelement(
                labels_g, '%d' % circle, ('x',0), ('y',radius-3))
            #labels_g.appendChild(t)

    def starsize(self, hipid):
        """Given a Hipparcos ID, return the radius of the size.

        This will return, for a given star, what the size of the
        circle of the drawing of said star should be.  If there is no
        star with that ID, then 0 is returned."""
        #if hipid<0 or len(self.hip_stars)<=hipid: return 0
        s = self.hip_stars[hipid]
        if s==None: return 0
        #return self.zerosize*(.8**(s[1]))
        #return self.zerosize-s[1]-2
        return self.dimmest_mag-s[1]+1

    def visible(self, hipid):
        """Given a Hipparcos ID, return if it is visible in this chart."""
        s = self.hip_stars[hipid]
        if s[3]<min(self.inner_dec, self.outer_dec): return False
        return s[3]<=max(self.inner_dec, self.outer_dec)

    def stars(self, magnitude=20):
        """Draw stars up to a certain magnitude."""
        # Get the stars that are visible within this chart.
        thestars = []
        for s in self.hip_stars:
            if not s: continue
            hip_id, mag, ra, dec, bv = s
            if mag>magnitude: continue
            if dec<min(self.inner_dec, self.outer_dec): continue
            if dec>max(self.inner_dec, self.outer_dec): continue
            thestars.append(s)
        # This should sort them by increasing magnitude (brightest first).
        thestars.sort(key=lambda a:a[1])
        if not thestars: return
        # Set the least bright magnitude.
        self.dimmest_mag = math.floor(thestars[-1][1])
        # Create the star group.
        star_g = self.make_element(self.centered, 'g', (
            'stroke', 'none'), ('fill', 'black'), (
            'clip-path', 'url(#innerClipPath)'))
        for hip_id, mag, ra, dec, bv in thestars:
            x, y = self.radec2xy(ra, dec)
            self.make_element(star_g, 'circle', (
                'cx', x), ('cy', y), ('r', self.starsize(hip_id)))

    def constellations(self):
        """Draw in the constellations."""
        # Create the constellation group.
        constellation_g = self.make_element(self.centered, 'g', (
            'stroke-width', '.5'), ('stroke-dasharray', '3,3'), (
            'clip-path', 'url(#innerClipPath)'))
        padding = 2
        for name, hips in constellation.constellations:
            for hip1, hip2 in hips:
                star1, star2 = self.hip_stars[hip1], self.hip_stars[hip2]
                if self.visible(hip1) or self.visible(hip2):
                    r1, r2 = self.starsize(hip1), self.starsize(hip2)
                    x1, y1 = self.radec2xy(star1[2], star1[3])
                    x2, y2 = self.radec2xy(star2[2], star2[3])
                    dx, dy = x2-x1, y2-y1
                    dd = math.sqrt(dx*dx+dy*dy)
                    if dd: dx, dy = dx/dd, dy/dd
                    x1 += dx*(r1+padding)
                    y1 += dy*(r1+padding)
                    x2 -= dx*(r2+padding)
                    y2 -= dy*(r2+padding)
                    line = self.make_element(constellation_g, 'line', (
                        'x1',x1), ('y1',y1), ('x2',x2), ('y2',y2))

    def __str__(self):
        return self.svg.toprettyxml()

if __name__ == "__main__":
    loc = body.ithaca
    time = body.Convert.date2day_number()
    nightb, nighte = body.sun.night(time)
    print time, nightb, nighte
    sr_nightb, sr_nighte = loc.sidereal(nightb), loc.sidereal(nighte)
    # Make the mini chart.
    chart = Chart(loc, sr_nightb, sr_nighte)
    chart.stars(5)
    chart.constellations()
    f = file('test.svg', 'w')
    print >> f, chart
    f.close()
