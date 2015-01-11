"""

Typically, all arguments of time are in terms of day numbers, a float
indicating the number of solar days since the start of the J2000
epoch.  For convenience sake, a convenience function
Convert.date2day_number is provided.

All angular arguments are typically in units of A rather irritating
situation in astronomy is inconsistency in

Locations are specified through the use of the 

This is code to compute the positions of the planets.

How to compute planetary positions
Paul Schlyter, Stockholm, Sweden, pausch@stjarnkimlen.se
http://stjarnhimlen.se/comp/ppcomp.html"""

import datetime, math, pdb

# When did the J2000 epoch start precisely, in UTC time?
j2000_zero_date = datetime.datetime(2000, 1, 1, 11, 58, 55, 816000)

class Convert:
    """Conversions relevant to astronomical calculation."""
    @classmethod
    def fractional(self, hours, minutes, seconds):
        "Fractional hours (degs) given hours (degrees), minutes, and seconds."
        return (seconds+60*(minutes+60*hours))/3600.0

    @classmethod
    def date2day_number(self, target_date = None):
        """Return the day number for a given datetime.datetime object.

        The day number is the number of days from the J2000 epoch.  If
        no target date is given, then assume we want the date from now."""
        # If no date is specified, make it right now.
        if target_date == None: target_date = datetime.datetime.utcnow()
        # Calculate and return the day number.
        y, m, d = target_date.year, target_date.month, target_date.day
        ut = target_date.hour+target_date.minute/60.0+target_date.second/3600.0
        dn = 367*y - 7*(y+(m+9)/12)/4 + 275*m/9 + d - 730530 + ut/24.0
        return dn

    @classmethod
    def day_number2date(self, d):
        """Given a day number, get the corresponding UTC datetime.datetime."""
        if d==None: return None
        zero_date = datetime.datetime(1999, 12, 31, 0, 0, 0)
        return zero_date + datetime.timedelta(d)

    @classmethod
    def sidereal2day(self, s):
        """For a radian sidereal interval, compute the number of solar days.

        A sidereal interval is a measure of time, with 2*pi being
        equal to one sidereal day (which is itself about four minutes
        less than a solar day).  For a measure of time in a sidereal
        interval, return the equivalent number of solar days."""
        return s/6.3003880986133822

class Format:
    """Formattings relevant to astronomical display."""
    @classmethod
    def rad2h(self, r):
        """An hour/minute/second string display of a given radian angle."""
        h = (math.degrees(r)%360.0)/15
        return '%2d:%02d:%04.1f' % (int(h), int(h*60)%60, (h*3600)%60)

    @classmethod
    def rad2d(self, r):
        """A degree/minute/second string display of a given radian angle."""
        d = ((math.degrees(r)+90)%360.0)-90
        ad = abs(d)
        return "%s%3d %02d'%04.1f\"" % (
            '+-'[d<0],int(ad), int(ad*60)%60, (ad*3600)%60)

class Location:
    """A location on Earth, in longitude and latitude."""

    # This holds the last location declared.
    last = None
    
    def __init__(self, longitude, latitude):
        """Create a location with a longitude and latitude, in radians."""
        self.lo, self.la = float(longitude), float(latitude)
        Location.last = self

    def sidereal(self, d=None):
        """Compute the local sidereal time at a target day number.

        Though we measure traditionally days in terms of periodicity
        of the suns apparent movement (the solar day), because the
        earth is moving around the sun, a single rotation of the earth
        is not identical to the period of the sun's movement.  As a
        result, stars appear to rise and set at different times of the
        day.  The notion of a 'sidereal' day is in terms of the
        periodicity of the earth's actual rotation, and consequently
        of the stars movement.

        This function, given a day number (as returned by
        Convert.date2day_number, perhaps) shall give the local (i.e.,
        longtidue adjusted) sidereal time at this location.

        This is relevant to astronomical calculations because when an
        object's right ascension is equal to the sidereal time at a
        location, it is at its zenith (or nadir) at that location.

        Note that the sidereal time is also, in the northern and
        southern hemisphere, the RA coordinate at due south or north,
        respectively."""
        if d==None: d=Body.default_day_number
        lst = math.radians(98.9818 + 0.985647345*d + (d%1)*360) + self.lo
        return lst % (2.0*math.pi)

    def equatorial(self, az, alt, d=None):
        """Return equatorial coordinates for azimuth and altitude coordinates.

        For a given location in the sky given the azimuth and altitude
        coordinates, this will return the corresponding equatorial
        coordinates.  Azimuth are the radians from north (moving
        initially to the east).  The returned values are (ra, dec)."""
        if d==None: d=Body.default_day_number
        # Method of conversion taken from
        # http://en.wikipedia.org/wiki/Horizontal_coordinate_system
        sindec = math.sin(self.la)*math.sin(alt)+\
                 math.cos(self.la)*math.cos(alt)*math.cos(az)
        cosdec_cosha = math.cos(self.la)*math.sin(alt)-\
                       math.sin(self.la)*math.cos(alt)*math.cos(az)
        cosdec_sinha = -math.sin(az)*math.cos(alt)
        ha = math.atan2(cosdec_sinha, cosdec_cosha)
        ra = self.sidereal(d) - ha
        dec = math.atan2(sindec, (cosdec_cosha**2+cosdec_sinha**2)**.5)
        return ra, dec

cf = Location(math.radians(-Convert.fractional(81,23,19)),
              math.radians(Convert.fractional(41,25,52)))

ithaca = Location(math.radians(-Convert.fractional(76,30,0)),
                  math.radians(Convert.fractional(42,26,36)))

redmond = Location(math.radians(-Convert.fractional(122,7,26)),
                   math.radians(Convert.fractional(47,40,10)))

alert = Location(math.radians(-Convert.fractional(62,20,20)),
                 math.radians(Convert.fractional(82,30,05)))

# The classes for representing stellar bodies.

class Body:
    """This is a sky element.
    
    Also included is the name, which allows one to look up
    previous elements by indexing on the name.

    The intension is for this to be subclassed.  The subclass should
    implement some of the methods that are here.  At minimum, an
    equatorial method that, given the day number, returns the
    equatorial coordinates, should be included."""
    # The default day number.  Initialized to be the current time, bu
    # can be reset.
    default_day_number = Convert.date2day_number()
    # Map of names to existing orbital elements already declared.
    existing, existing_names = {}, []

    def __init__(self,names):
        if isinstance(names, str): names = [names]
        self._names = names
        if False and name:
            if name in self.existing:
                self.existing_names.remove(name)
            self.existing[name] = self
            self.existing_names.append(name)

    def describe(self):
        return ' / '.join(self._names)

    def names(self):
        return self._names

    def equatorialt(self, d=None, loc=None):
        """Equatorial coordinates as viewed on an earthbound location.

        Geocentric equatorial coordinates give an idea of an object's
        location in the sky, but closeby objects may be subject to
        parallax relative to an observer's location on the.  These
        'surface' coordinates are instead called the topocentric
        equatorial coordinates.

        The default behavior of this is to just return the regular
        geocentric equatorial coordiantes, which is certainly
        acceptable for anything beyond the solar system, and generally
        considered acceptable for anything other than the moon.

        Subclasses may override this method to factor in any
        appropriate desired adjustments."""
        return self.equatorial(d)

    def azalt(self, d=None, loc=None):
        """Return the azimuth and altitude in radians.

        This will extract the object's equatorial coordinates and, for
        the given day number and location, return the alt-az of that
        object."""
        if loc==None: loc=Location.last
        ra, dec = self.equatorialt(d, loc)
        # Get local sidereal time, and compute the hour angle.
        ha = loc.sidereal(d) - ra
        # Finally, compute the altitude and azimuth.
        alt = math.asin( math.sin(loc.la)*math.sin(dec)
                        +math.cos(loc.la)*math.cos(dec)*math.cos(ha))
        az = math.atan2(-math.sin(ha), -math.cos(ha)*math.sin(loc.la)
                        +math.tan(dec)*math.cos(loc.la))
        return az%(2.0*math.pi), (alt+math.pi)%(2.0*math.pi)-math.pi

    def meridian(self, offset=0, d=None, loc=None):
        """Return the day number of a meridian of the object.

        If offset=0, then this will find either the last or next
        meridian, with a tendency towards the closest meridian.  If
        offset=-1 or 1, then it will find the last or next meridian
        respectively.  Offset can also be <-1 or >1, in which case it
        will attempt to look to futher offsets.  If the offset is
        positive, then it will attempt to find the

        Using high offsets will be more inaccurate if the equatorial
        coordinates of an object change over time.  It would be
        reasonably easy to make this much more if not perfectly
        accurate, but given that the current intention is that
        magnitude of the offset never exceed anything more than what
        one could count on one finger, for the interest in simplicty
        this inexactness is tolerated."""
        if loc==None: loc=Location.last
        if d==None: d=self.default_day_number
        nd = d
        while True:
            while True:
                lst = loc.sidereal(nd)
                ra, dec = self.equatorialt(nd, loc)
                tds = ((ra-lst)+math.pi)%(2*math.pi) - math.pi
                if abs(tds)<1e-8: break
                nd += Convert.sidereal2day(tds)
            trueoffset = -1 if nd < d else 1
            if offset==0 or trueoffset==offset: break
            # Hmmm, by what period should we change?
            period=offset-(max,min)[offset<0](0,trueoffset)
            nd += period*Convert.sidereal2day(2*math.pi)
            offset = 0 # Assume it'll be right after this.
        return nd

    def nadir(self, offset=0, d=None, loc=None):
        """Return the day number the object is closest to nadir.

        If offset=0, then this will find either the last or next
        meridian, with a tendency towards the closest meridian.  If
        offset=-1 or 1, then it will find the last or next meridian
        respectively.  Offset can also be <-1 or >1, in which case it
        will attempt to look to futher offsets.  If the offset is
        positive, then it will attempt to find the

        Using high offsets will be more inaccurate if the equatorial
        coordinates of an object change over time.  It would be
        reasonably easy to make this much more if not perfectly
        accurate, but given that the current intention is that
        magnitude of the offset never exceed anything more than what
        one could count on one finger, for the interest in simplicty
        this inexactness is tolerated."""
        if loc==None: loc=Location.last
        if d==None: d=self.default_day_number
        nd = self.meridian(offset, d, loc)
        while True:
            lst = loc.sidereal(nd)
            ra, dec = self.equatorialt(nd, loc)
            tds = (ra-lst+15*math.pi/16)%(2*math.pi) - 31*math.pi/16
            if abs(tds)<1e-8: break
            nd += Convert.sidereal2day(tds)
        return nd

    def rise(self, h=0, offset=0, d=None, loc=None):
        """Return the day number of the time the object rises.

        The h parameter is the angle above the horizon which
        constitutes a rise.  Which is appropriate will depend on the
        circumstances.

        The offset argument indicates that we should find the rise
        directly before the meridian that would be found by feeding
        this offset parameter into the function.  If a rise associated
        with the meridian cannot be detected for the object (e.g.,
        because the object is always above or below the horizon), this
        function shall return None."""
        if loc==None: loc=Location.last
        if d==None: d=self.default_day_number
        h = math.radians(h)
        # Get the day number of the particular meridian of this object we want.
        dz = self.meridian(offset, d, loc)
        risetime = nd = dz
        tds=0
        while True:
            lst = loc.sidereal(nd)
            ra, dec = self.equatorialt(nd, loc)
            # Get the apparently closest meridian time's sidereal
            # offset.  If the object is stationary in the sky, we want
            # one which is after or on this time, given that we are
            # searching for a rise.  On the other hand, in iteration,
            # with a object rapidly changing in the sky, if its .  So,
            # we encourage postivity far more than we would by
            # considering offsets in the range [-pi/16, 31pi/16],
            # instead of the old [-pi,pi].
            tds = ((ra-lst)+math.pi/16)%(2*math.pi) - math.pi/16
            # Calculate the sidereal offset of the time.
            cLHA = math.sin(h)-math.sin(loc.la)*math.sin(dec)
            cLHA /= math.cos(loc.la)*math.cos(dec)
            if -1>cLHA or cLHA>1:
                return None # Never crosses this.
            # Now... we must adjust the new day.
            LHA = math.acos(cLHA)
            nd += Convert.sidereal2day(tds-LHA)
            if abs(nd-risetime)<1.0/86400: break
            risetime = nd
        return risetime

    def set(self, h=0, offset=0, d=None, loc=None):
        """Return the day number of the time the object sets.

        This function behaves similarly to the rise function, except
        that offset indicates which meridian that we should look for a
        set after, not before as in the rise function."""
        if loc==None: loc=Location.last
        if d==None: d=self.default_day_number
        h = math.radians(h)
        # Get the day number of the particular meridian of this object we want.
        dz = self.meridian(offset, d, loc)
        risetime = nd = dz
        tds=0
        while True:
            lst = loc.sidereal(nd)
            ra, dec = self.equatorialt(nd, loc)
            # Get the apparently relevant meridian time's sidereal
            # offset.  For similar reasons to the rise calculation, we
            # consider meridian offsets [-31pi/16,pi/16].
            tds = ((ra-lst)+31*math.pi/16)%(2*math.pi) - 31*math.pi/16
            # Calculate the sidereal offset of the time.
            cLHA = math.sin(h)-math.sin(loc.la)*math.sin(dec)
            cLHA /= math.cos(loc.la)*math.cos(dec)
            if -1>cLHA or cLHA>1:
                return None # Never crosses this.
            # Now... we must adjust the new day.
            LHA = math.acos(cLHA)
            nd += Convert.sidereal2day(tds+LHA)
            if abs(nd-risetime)<1.0/86400: break
            risetime = nd
            #print risetime; import time; time.sleep(.1)
        return risetime

    def events(self, s, e, h=None, loc=None):
        """This method will return a list of events between two times.

        The two times are a range start 's' and end 'e', specified as
        usual in day numbers.  'h' is how high in the sky the object
        should be to be considered up.  Without any information, it is
        half of the object's apparent diameter at the start time.

        Given a start and an end time, this will return a list of all
        relevant rises, meridians, and sets of the object.  The list
        contains one item for each event, where the event is a two
        element tuple containing, first, the day number of the event,
        and then the single character code for the event.  Codes are:

        N    - The object is closest to nadir.
        R    - The object rose.
        M    - The object reached meridian.
        S    - The object set.
        s, e - Start and end times within a period the object is up.

        Note that if the start time occurs before a set, or the end
        time occurs after a rise, that the rise and meridian
        associated with that set (or merdian and set associated with
        that rise) shall be included in the event list, even if they
        are outside the interval.  In this 'intercepting' case, there
        shall also be added events with codes 's' and 'e' (note the
        lower case) corresponding to the start and end of the defined
        interval.  These sort of boundary events are added only when
        the object was up during the start."""
        assert s < e
        useu = False
        if h==None: h=self.apparent_diameter(s)/2.0
        events = []
        # Find the first meridian before the start time.
        events.extend([
            (self.meridian(-1, s, loc),'M'), (self.set(h, -1, s, loc),'S'),
            (self.nadir(1, s, loc),'N'), (self.rise(h, 1, s, loc),'R'),
            (self.meridian(1, s, loc),'M'), (self.set(h, 1, s, loc),'S')])
        m = events[-2][0]
        while m < e:
            m+=.1
            events.extend([
                (self.nadir(1, m, loc),'N'), (self.rise(h, 1, m, loc),'R'),
                (self.meridian(1, m, loc),'M'), (self.set(h, 1, m, loc),'S')])
            m = events[-2][0]
        events = [(dn,c) for dn,c in events if dn!=None]
        events.append((s,'s'))
        events.append((e,'e'))
        events.sort()
        # Now, run through the event list.
        nevents = []
        inday = self.azalt(events[0][0], loc)>=h
        ending = False
        for daynum, code in events:
            if code=='s':
                if inday:
                    del nevents[:-1]
                    nevents.append((daynum,code))
                else:
                    nevents = []
            elif code=='e':
                if not inday: break
                ending = True
                nevents.append((daynum,code))
            elif code in 'MNs' and inday:
                nevents.append((daynum,code))
                if ending: break
            if code in 'RS':
                nevents.append((daynum,code))
                if ending: break
                inday = code=='R'
        return nevents

    def apparent_diameter(self, d=None):
        """Radian arc of the object's diameter from a geocentric observer.

        Defaults to 0 to indicate a pinpoint object of no discernable
        radius, as in a star.  Subclasses may override for more
        intelligent behavior, as one might like for a solar system
        object or a deep-sky object."""
        return 0.0

    @classmethod
    def obliquity(cls, d=None):
        """Ecliptic obliquity for the earth's orbit at a target day number."""
        if d==None: d=cls.default_day_number
        return math.radians(23.4393 - 3.563E-7*d)

class DeepSkyBody(Body):
    """This contains a deep sky object's key parameters."""
    def __init__(self, name, mag, ra, dec, pmra=0.0, pmdec=0.0, diam=0.0):
        """Creates a deep sky body.

        Arguments indicate the name, a given magnitude, right
        ascension in hours and declination in degrees, proper motion
        of the right ascension and declination per year in
        milli-arc-seconds, and apparent diameter in the sky in
        degrees."""
        Body.__init__(self,name)
        # Make the right ascension and declination.
        self._mag, self._diam = mag, diam
        self._ra, self._dec = math.radians(ra*15), math.radians(dec)
        f = 36e6 * 365.25 # marcsec/year translated into degrees/day
        self._pmra, self._pmdec = math.radians(pmra/f), math.radians(pmdec/f)
        self._diam = math.radians(diam)

    def equatorial(self, d=None):
        """Return the radian RA and Dec coordinates of the object."""
        if d==None: d=self.default_day_number
        return self._ra+self._pmra*d, self._dec+self._pmdec*d

    def apparent_diameter(self, d=None):
        """Get the radian apparent diameter of the object from earth's POV."""
        return self._diam

    def magnitude(self, d=None):
        """Get the apparent magnitude of the object."""
        return self._mag

class OrbitalBody(Body):
    """This contains a solar system's body's characteristic.

    Note that returned angles are in radians rather than degrees.
    This may be confusing, because angles in the initialization
    functions are given in degrees, but this is for convenience in
    calculation and entry, respectively."""

    def __init__(self,name,dequ,dpol,mprime,mphase,
                 N,i,w,a,e,M,Nd=0.0,id=0.0,wd=0.0,ad=0.0,ed=0.0,Md=0.0):
        """Initializes a structure containing a body's characteristics.

        For apparent diameter, dequ and dpol are the apparent diameter
        in arcseconds from 1 distance unit away (usually AU), both at
        the equator and the pole respectively.  (The equatorial
        diameter is presumably no smaller than the polar diameter.)

        For magnitude, mprime and mphase are the constant and
        phase-multiplied terms of the magnitude equation.

        There are also orbital elements one enters for a body.  The
        primary orbital elements are,
        N - longitude of the ascending node
        i - inclination of the ecliptic
        w - argument of perihelion
        a - semi-major axis, or mean distance from sun (in AU)
        e - eccentricity (0=circle, 0-1=ellipse, 1=parabola)
        M - mean anomaly (0 at perhelion, increases uniformly with time)

        There are also associated Nd, id, etc., arguments for each of
        these elements.  Why?  Many of these quantities differ across
        days, so when retrieving the elements through the 'elements'
        method , one must specify a day number.  For example, the
        elements w, e, etc., when returned will be w+wd*d, e+ed*d,
        etc.  All of these elements have optional 'day progression'
        quantities (defaulting to 0.0) which allows one to capture
        this change across time."""
        Body.__init__(self,name)
        # Input apparent diameter elements.
        self._dequ, self._dpol = (math.radians(i/3600.0) for i in (dequ, dpol))
        # Input apparent magnitude elements.
        self._mprime, self._mphase = mprime, mphase
        # Input orbital elements.
        self._N,self._Nd,self._i,self._id = N, Nd, i, id
        self._w,self._wd,self._a,self._ad = w, wd, a, ad
        self._e,self._ed,self._M,self._Md = e, ed, M, Md

    def name(self):
        """Returns the name of this object."""
        return self._name

    def N(self, d=None):
        """Longitude of the ascending node."""
        if d==None: d=self.default_day_number
        return math.radians((self._N+self._Nd*d)%360.0)

    def i(self, d=None):
        """Inclination of the ecliptic."""
        if d==None: d=self.default_day_number
        return math.radians((self._i+self._id*d)%360.0)

    def w(self, d=None):
        """Argument of the perihelion."""
        if d==None: d=self.default_day_number
        return math.radians((self._w+self._wd*d)%360.0)

    def a(self, d=None):
        """The semi-major axis, or mean distance from the sun."""
        if d==None: d=self.default_day_number
        return self._a+self._ad*d

    def e(self, d=None):
        """The eccentricity (0=circle, 0-1=ellipse, 1=parabola)."""
        if d==None: d=self.default_day_number
        return self._e+self._ed*d

    def M(self, d=None):
        """Longitude of the ascending node."""
        if d==None: d=self.default_day_number
        return math.radians((self._M+self._Md*d)%360.0)

    def w1(self, d=None):
        """Longitude of the perihelion."""
        if d==None: d=self.default_day_number
        return math.radians((self._N+self._w+(self._Nd+self._wd)*d)%360.0)

    def L(self, d=None):
        """Mean longitude."""
        if d==None: d=self.default_day_number
        return math.radians((self._M+self._N+self._w+(
            self._Md+self._Nd+self._wd)*d)%360.0)

    def q(self, d=None):
        """Perihelion distance."""
        if d==None: d=self.default_day_number
        return self.a(d)*(1.0-self.e(d))

    def Q(self, d=None):
        """Aphelion distance."""
        if d==None: d=self.default_day_number
        return self.a(d)*(1.0+self.e(d))

    def P(self, d=None):
        """Orbital period in years."""
        if d==None: d=self.default_day_number
        return self.a**1.5

    def T(self, d=None):
        """Epoch of M."""
        if d==None: d=self.default_day_number
        return 2000.0 # Um, just nothing.

    def V(self, d=None):
        """True anomaly (angle between position and perihelion)."""
        if d==None: d=self.default_day_number
        e, E = self.e(d), self.E(d)
        V = 2.0 * math.atan(math.sqrt((1+e)/(1-e))*math.tan(0.5*E))
        return V%(2*math.pi)

    def E(self, d=None, tol=None):
        """Eccentric anomaly, calculated to a tolerance."""
        if d==None: d=self.default_day_number
        if tol==None: tol=1e-6
        M, e = self.M(d), self.e(d)
        En = M+e*math.sin(M)*(1.0+e*math.cos(M))
        while True:
            E = En
            En = E - (E-e*math.sin(E)-M)/(1.0-e*math.cos(E))
            if abs(En-E) <= tol: break
        E = En%(2*math.pi)
        return E

    def R(self, d=None):
        """Get the helocentric radius of the object."""
        if d==None: d=self.default_day_number
        a, e, V = self.a(d), self.e(d), self.V(d)
        R = a*(1.0-e*e) / (1.0+e*math.cos(V))
        return R

    def lon_lat_correction(self, d=None):
        """Return longitude/latitude perturbation offset.  Default is 0,0.

        When calculating the ecliptic longitude and latitude, there
        are perturbations that can arise due to the influence of
        bodies other than the sun.  This will differ from object to
        object.

        This function is used in computing the ecliptic longitude and
        latitude, which applys this correction for perturbation.  The
        intention is the peculiar objects with strange perturbations
        will produce subclasses that override this function."""
        return 0.0, 0.0

    def ecliptic_lonlat(self, d=None):
        """Return the ecliptic longitude and latitude."""
        # First calculate pure mathematical unperturbed heliocentric coords.
        N, V, w, i, r = self.N(d), self.V(d), self.w(d), self.i(d), self.R(d)
        xh=r*(math.cos(N)*math.cos(V+w)-math.sin(N)*math.sin(V+w)*math.cos(i))
        yh=r*(math.sin(N)*math.cos(V+w)+math.cos(N)*math.sin(V+w)*math.cos(i))
        zh=r*(math.sin(V+w)*math.sin(i))
        # Calculate the longitude, latitude of the ecliptic (from vernal pt).
        lonecl = math.atan2(yh, xh)
        latecl = math.atan2(zh, math.hypot(xh, yh))
        # Apply the perturbation corrections to the ecliptic coordinates.
        lonc, latc = self.lon_lat_correction(d)
        return lonecl+lonc, latecl+latc

    def xyzh(self, d=None):
        """Return the heliocentric (x,y,z) coordinates."""
        r, (lonecl, latecl) = self.R(d), self.ecliptic_lonlat(d)
        xh = r*math.cos(lonecl)*math.cos(latecl)
        yh = r*math.sin(lonecl)*math.cos(latecl)
        zh = r*math.sin(latecl)
        return xh, yh, zh

    def earth_distance(self, d=None):
        """Return the earth's distance from this object."""
        # Otherwise we must offset properly.
        xh, yh, zh = self.xyzh(d)
        xs, ys = self.sun.xys(d)
        return math.sqrt((xh+xs)*(xh+xs)+(yh+ys)*(yh+ys)+zh*zh)

    def parallax(self, d=None):
        """Return the apparent size in radians of earth's equatorial radius.

        Note that this is different from the more common meaning of
        parallax, which is applied to the shift one sees in an entire
        half-orbital period rather than a single half-orbital
        rotation."""
        dist = self.earth_distance(d)
        if isinstance (self, OrbitalBodyMoon):
            return math.asin( 1.0 / dist )
        else:
            # The radius of the earth is about 4.2487e-5 AUs.
            return math.asin( 4.2487e-5 / dist )

    def equatorial(self, d=None):
        """Return the geocentric equatorial coordinates (RA and Dec).

        As consistent, the returned floats are in radians."""
        # Get the heliocentric coordinates of this object.
        xh, yh, zh = self.xyzh(d)
        if not isinstance(self, OrbitalBodyMoon):
            # Get the x_s, y_s solar position.
            xs, ys = self.sun.xys(d)
            # Now make it geocentric through a simple offset.
            xg, yg, zg = xh+xs, yh+ys, zh
        else:
            xg, yg, zg = xh, yh, zh
        # Now rotate into the equatorial plane.
        ecl = self.obliquity(d)
        xe, ye, ze = xg, yg*math.cos(ecl)-zg*math.sin(ecl), \
                     yg*math.sin(ecl)+zg*math.cos(ecl)
        # Now, return the right ascension and declination.
        ra, dec = math.atan2(ye, xe), math.atan2(ze, math.hypot(xe, ye))
        return ra, dec

    def equatorialt(self, d=None, loc=None):
        """Return the topocentric (surface) equatorial coordinates.

        This method is desirable because for some very close objects,
        like the moon, the difference between different locations on
        the earth will vary quite a bit.

        This method may take not only a time, but a location.  If
        location is left to its default of None, then the last
        declared Location object, Location.last, is used instead."""
        if loc==None: loc=Location.last
        # Get the geocentric equatorial things.
        ra, dec = self.equatorial(d)
        # Get the local sidereal time.
        lst = loc.sidereal()
        # What is the parallax of the object, i.e., apparent size of
        # the equatorial radius of earth from our distant object?
        par = self.parallax(d)
        # Compute geocentric latitude, and distance rho from the earth center.
        gclat = loc.la-0.003358*math.sin(2*loc.la)
        rho = 0.99833+0.00167*math.cos(2*loc.la)
        # Compute the hour angle of the object in question.
        ha = lst - ra
        # Get the auxiliary angle g.
        g = math.atan( math.tan(gclat) / math.cos(ha) )
        # Now, make the topocentric RA and Dec.
        cdec = math.cos(dec)
        if cdec > 1e-3:
            topRa = ra - par*rho*math.cos(gclat)*math.sin(ha)/cdec
        else:
            topRa = ra # No correction, formula invalid anyway...
        if abs(g) > 1e-4:
            topDec = dec - par*rho*math.sin(gclat)*math.sin(g-dec)/math.sin(g)
        else:
            topDec = dec - par*rho*math.sin(g-dec)*math.cos(ha)
        return topRa, topDec

    # Methods about the appearance of the planet.

    def elongation(self, d=None):
        """Get elongation, planet-sun angular distance from earth's POV."""
        s, r, R = self.sun.R(d), self.R(d), self.earth_distance(d)
        return math.acos((s*s+R*R-r*r)/(2*s*R))

    def apparent_diameter(self, d=None):
        """Get the apparent diameter of the object from earth's POV."""
        return (self._dequ + self._dpol) / (2*self.earth_distance(d))

    def phase_angle(self, d=None):
        """Return the phase angle."""
        s, r, R = self.sun.R(d), self.R(d), self.earth_distance(d)
        return math.acos((r*r+R*R-s*s) / (2*r*R))

    def phase(self, d=None):
        """Return the phase."""
        return (1.0-math.cos(self.phase_angle(d)))/2.0

    def magnitude(self, d=None):
        """Return the apparent magnitude of the object."""
        r, R = self.R(d), self.earth_distance(d)
        fv = math.degrees(self.phase_angle(d))
        return self._mprime + 5.0*math.log10(r*R) + self._mphase*fv

class OrbitalBodySun(OrbitalBody):
    def __init__(self, *args):
        OrbitalBody.__init__(self, *args)
        OrbitalBody.sun = self

    def xys(self, d=None):
        e, E, R = self.e(d), self.E(d), self.R(d)
        xy, yv = math.cos(E)-e, math.sqrt(1.0-e*e)*math.sin(E)
        lonsun = self.V(d) + self.w(d)
        xs, ys = R*math.cos(lonsun), R*math.sin(lonsun)
        return xs, ys

    def earth_distance(self, d=None):
        return self.R(d)

    def equatorial(self, d=None):
        (xs, ys), ecl = self.xys(d), self.obliquity(d)
        xe, ye, ze = xs, ys*math.cos(ecl), ys*math.sin(ecl)
        ra, dec = math.atan2(ye, xe), math.atan2(ze, math.hypot(xe, ye))
        return ra, dec

    def elongation(self, d=None):
        """Elongation for the sun is obviously always 0."""
        return 0.0

    def phase(self, d=None): return 0.0 # Sun doesn't have 'phases' as such...

    def magnitude(self, d=None):
        # Use the simple distance formula for visual magnitude.
        return 4.83 + 5*(math.log10(self.R(d)/206264.806)-1)

    def night(self, d=None, loc=None, h=-15.0):
        """Returns day numbers of the next (or current) night.

        This is a convenience method for the sun that returns the
        two-element tuple of the next or current night from the day
        number indicated from the current d.  If d is before a sunrise
        and after a sunset (at night), or else d is before some sunset
        but before the next sunrise (in day), then return that
        (sunset, sunrise) pair.

        One may define 'night' in different ways by the h argument.
        The default is for h=-15, so that night starts at amateur
        astronomical twilight.  One may also define h=-6 for civil
        twilight."""
        if d==None: d=self.default_day_number
        # Get the rise time of the next solar meridian.
        risedn = self.rise(h, 1, d, loc)
        if risedn < d:
            # If the rise is before now, then we should return the next
            # sunset, and the sunrise after that.
            return self.set(h, 1, d, loc), self.rise(h, 2, d, loc)
        # If the rise if after now, then we should return the last
        # sunset to that next sunrise.
        return self.set(h, -1, d, loc), risedn
        

class OrbitalBodyMoon(OrbitalBody):
    def earth_distance(self, d=None):
        """Return the earth's distance from this object."""
        sun = self.sun
        Ms,Mm,Nm,ws,wm = sun.M(d), self.M(d), self.N(d), sun.w(d), self.w(d)
        Ls, Lm = Ms+ws, Mm+wm+Nm
        D = Lm - Ls
        return self.R(d)-0.58*math.cos(Mm-2*D)-0.46*math.cos(2*D)

    def lon_lat_correction(self, d=None):
        sun = self.sun
        Ms,Mm,Nm,ws,wm = sun.M(d), self.M(d), self.N(d), sun.w(d), self.w(d)
        Ls, Lm = Ms+ws, Mm+wm+Nm
        D = Lm - Ls
        F = Lm - Nm
        londelta = -1.274 * math.sin(Mm - 2*D) \
                   +0.658 * math.sin(2*D) \
                   -0.186 * math.sin(Ms) \
                   -0.059 * math.sin(2*Mm - 2*D) \
                   -0.057 * math.sin(Mm - 2*D + Ms) \
                   +0.053 * math.sin(Mm + 2*D) \
                   +0.046 * math.sin(2*D - Ms) \
                   +0.041 * math.sin(Mm - Ms) \
                   -0.035 * math.sin(D) \
                   -0.031 * math.sin(Mm + Ms) \
                   -0.015 * math.sin(2*F - 2*D) \
                   +0.011 * math.sin(Mm - 4*D)
        latdelta = -0.173 * math.sin(F - 2*D) \
                   -0.055 * math.sin(Mm - F - 2*D) \
                   -0.046 * math.sin(Mm + F - 2*D) \
                   +0.033 * math.sin(F + 2*D) \
                   +0.017 * math.sin(2*Mm + F)
        return math.radians(londelta), math.radians(latdelta)

    def elongation(self, d=None):
        (slon, slat), (mlon, mlat) = (
            b.ecliptic_lonlat(d) for b in (self.sun, self))
        return math.acos(math.cos(slon-mlon)*math.cos(mlat))

    def phase_angle(self, d=None):
        return math.pi - self.elongation(d)

    def magnitude(self, d=None):
        r, R = self.sun.R(d), self.earth_distance(d)/23450.0
        fv = math.degrees(self.phase_angle(d))
        # Use the cosine rule to make r into the true distance...
        oldr = r
        r = R*math.cos(fv)+math.sqrt(r*r-R*math.sin(fv)**2)
        return self._mprime+5.0*math.log10(r*R)+self._mphase*fv+4e-9*(fv**4)

class OrbitalBodyMercury(OrbitalBody):
    def magnitude(self, d=None):
        m,fv = OrbitalBody.magnitude(self,d), math.degrees(self.phase_angle(d))
        return m + 2.2e-13 * fv**6

class OrbitalBodyVenus(OrbitalBody):
    def magnitude(self, d=None):
        m,fv = OrbitalBody.magnitude(self,d), math.degrees(self.phase_angle(d))
        return m + 4.2e-7 * fv**3

class OrbitalBodyJupiter(OrbitalBody):
    def __init__(self, *args):
        OrbitalBody.__init__(self, *args)
        OrbitalBody.jupiter = self
    def lon_lat_correction(self, d=None):
        Mj, Ms = self.M(d), self.saturn.M(d)
        londelta = -0.332 * math.sin(2*Mj-5*Ms-1.1816) \
                   -0.056 * math.sin(2*Mj - 2*Ms + 0.36652) \
                   +0.042 * math.sin(3*Mj - 5*Ms + 0.36652) \
                   -0.036 * math.sin(Mj - 2*Ms) \
                   +0.022 * math.cos(Mj - Ms) \
                   +0.023 * math.sin(2*Mj - 3*Ms + 0.90757) \
                   -0.016 * math.sin(Mj - 5*Ms - 1.2043)
        return math.radians(londelta), 0.0

class OrbitalBodySaturn(OrbitalBody):
    def __init__(self, *args):
        OrbitalBody.__init__(self, *args)
        OrbitalBody.saturn = self
    def lon_lat_correction(self, d=None):
        Ms, Mj = self.M(d), self.jupiter.M(d)
        londelta = +0.812 * math.sin(2*Mj - 5*Ms - 1.1816) \
                   -0.229 * math.cos(2*Mj - 4*Ms - 0.03591) \
                   +0.119 * math.sin(Mj - 2*Ms - 0.05234) \
                   +0.046 * math.sin(2*Mj - 6*Ms - 1.2043) \
                   +0.014 * math.sin(Mj - 3*Ms + 0.55851)
        latdelta = -0.020 * math.cos(2*Mj - 4*Ms - 0.03591) \
                   +0.018 * math.sin(2*Mj - 6*Ms - 0.85521)
        return math.radians(londelta), math.radians(latdelta)
    def magnitude(self, d=None):
        if d==None: d=self.default_day_number
        m = OrbitalBody.magnitude(self,d)
        # First get the geocentric ecliptic longitude and latitude.
        (xs,ys),(xh,yh,zh) = self.sun.xys(d), self.xyzh(d)
        los = math.atan2(yh+ys, xh+xs)
        las = math.atan2(zh, math.hypot(xh+xs, yh+ys))
        ir, Nr = 0.48974, 2.9585+6.667e-7*d
        sB = math.sin(las)*math.cos(ir)-\
             math.cos(las)*math.sin(ir)*math.sin(los-Nr)
        return m + -2.6*abs(sB) + 1.2*sB**2

class OrbitalBodyUranus(OrbitalBody):
    def __init__(self, *args):
        OrbitalBody.__init__(self, *args)
    def lon_lat_correction(self, d=None):
        Mu, Mj, Ms = self.M(d), self.jupiter.M(d), self.saturn.M(d)
        londelta = +0.040 * math.sin(Ms - 2*Mu + 0.10472) \
                   +0.035 * math.sin(Ms - 3*Mu + 0.57596) \
                   -0.015 * math.sin(Mj - Mu + 0.34907)
        return math.radians(londelta), 0.0

class OrbitalBodyPluto(OrbitalBody):
    def __init__(self, *args):
        OrbitalBody.__init__(self, *(args + ((0.0,)*12)))
    def P(self, d=None):
        if d==None: d=self.default_day_number
        return math.radians(238.95 + 0.003968789*d)
    def S(self, d=None):
        if d==None: d=self.default_day_number
        return math.radians(50.03 + 0.033459652*d)
    def R(self, d=None):
        if d==None: d=self.default_day_number
        P = self.P(d)
        r = 40.72 \
            + 6.68 * math.sin(P) + 6.90 * math.cos(P) \
            - 1.18 * math.sin(2*P) - 0.03 * math.cos(2*P) \
            + 0.15 * math.sin(3*P) - 0.14 * math.cos(3*P)
        return r
    def xyzh(self, d=None):
        if d==None: d=self.default_day_number
        P, S, r = self.P(d), self.S(d), self.R(d)
        lonecl = 238.9508 + 0.00400703 * d \
                 - 19.799 * math.sin(P) + 19.848 * math.cos(P) \
                 + 0.897 * math.sin(2*P) - 4.956 * math.cos(2*P) \
                 + 0.610 * math.sin(3*P) + 1.211 * math.cos(3*P) \
                 - 0.341 * math.sin(4*P) - 0.190 * math.cos(4*P) \
                 + 0.128 * math.sin(5*P) - 0.034 * math.cos(5*P) \
                 - 0.038 * math.sin(6*P) + 0.031 * math.cos(6*P) \
                 + 0.020 * math.sin(S-P) - 0.010 * math.cos(S-P)
        latecl = -3.9082 \
                 - 5.453 * math.sin(P) - 14.975 * math.cos(P) \
                 + 3.527 * math.sin(2*P) + 1.673 * math.cos(2*P) \
                 - 1.051 * math.sin(3*P) + 0.328 * math.cos(3*P) \
                 + 0.179 * math.sin(4*P) - 0.292 * math.cos(4*P) \
                 + 0.019 * math.sin(5*P) + 0.100 * math.cos(5*P) \
                 - 0.031 * math.sin(6*P) - 0.026 * math.cos(6*P) \
                 + 0.011 * math.cos(S-P)
        lonecl, latecl = math.radians(lonecl), math.radians(latecl)
        xh = r*math.cos(lonecl)*math.cos(latecl)
        yh = r*math.sin(lonecl)*math.cos(latecl)
        zh = r*math.sin(latecl)
        return xh, yh, zh

# Make a mapping of item names to orbital elements.
major_solar_system = []
sun = OrbitalBodySun(
    'Sun',1919.26,1919.26,-26.74,0.0,
    0.0, 0.0, 282.9404, 1.000000, 0.016709, 356.0470, 
    0.0, 0.0, 4.70935E-5, 0.0, -1.151E-9, 0.9856002585)
major_solar_system.append(sun)
major_solar_system.append(OrbitalBodyMoon(
    'Moon',1873.7*60,1873.7*60,0.23,0.026,
    125.1228, 5.1454, 318.0634, 60.2666, 0.054900, 115.3654, 
    -0.0529538083, 0.0, 0.1643573223, 0.0, 0.0, 13.0649929509))
major_solar_system.append(OrbitalBodyMercury(
    'Mercury',6.74,6.74,-0.36,0.026,
    48.3313, 7.0047, 29.1241, 0.387098, 0.205635, 168.6562, 
    3.24587E-5, 5.00E-8, 1.01444E-5, 0.0, 5.59E-10, 4.0923344368))
major_solar_system.append(OrbitalBodyVenus(
    'Venus',16.92,16.92,-4.34,0.013,
    76.6799, 3.3946, 54.8910, 0.723330, 0.006773, 48.0052, 
    2.46590E-5, 2.75E-8, 1.38374E-5, 0.0, -1.302E-9, 1.6021302244))
major_solar_system.append(OrbitalBody(
    'Mars',9.36,9.28,-1.51,0.016,
    49.5574, 1.8497, 286.5016, 1.523688, 0.093405, 18.6021, 
    2.11081E-5, -1.78E-8, 2.92961E-5, 0.0, 2.516E-9, 0.5240207766))
major_solar_system.append(OrbitalBodyJupiter(
    'Jupiter',196.94,185.08,-9.25,0.014,
    100.4542, 1.3030, 273.8777, 5.20256, 0.048498, 19.8950, 
    2.76854E-5, -1.557E-7, 1.64505E-5, 0.0, 4.469E-9, 0.0830853001))
major_solar_system.append(OrbitalBodySaturn(
    'Saturn',165.6,150.8,-9.0,0.044,
    113.6634, 2.4886, 339.3939, 9.55475, 0.055546, 316.9670, 
    2.38980E-5, -1.081E-7, 2.97661E-5, 0.0, -9.499E-9, 0.0334442282))
major_solar_system.append(OrbitalBodyUranus(
    'Uranus',65.8,62.1,-7.15,0.001,
    74.0005, 0.7733, 96.6612, 19.18171, 0.047318, 142.5905, 
    1.3978E-5, 1.9E-8, 3.0565E-5, -1.55E-8, 7.45E-9, 0.011725806))
major_solar_system.append(OrbitalBody(
    'Neptune',62.2,60.9,-6.90,0.001,
    131.7806, 1.7700, 272.8461, 30.05826, 0.008606, 260.2471,
    3.0173E-5, -2.55E-7, -6.027E-6, 3.313E-8, 2.15E-9, 0.005995147))
# This has relatively low error.
# The -1.0 figure for reflective magnitude of the planet came from
# Astronomy on the Personal Computer By Oliver Montenbruck, p 144, table 7.3
major_solar_system.append(OrbitalBodyPluto('Pluto', 3.1795, 3.1795, -1.0, 0.0))

# This, on the other hand, has high error.
#OrbitalBody(
#    'Pluto',
#    110.30347, 17.14175, 224.06676, 39.48168677, 0.24880766, 238.92881,
#    110.30347, 17.14175, 113.76329+180, 39.48168677, 0.24880766, 238.92881,
#    0.0, 0.0, 0.0, 0.0, 0.0, 360.0/90613.3055)

deep_sky_objects = []
deep_sky_objects.append(DeepSkyBody(
    'M31', 3.50, Convert.fractional(0,44,30), Convert.fractional(41,29,0),
    0, 0, Convert.fractional(0,1,17)))
deep_sky_objects.append(DeepSkyBody(
    'Polaris', 1.95,
    Convert.fractional(2,31,48.7), Convert.fractional(89,15,51.1)))
pol = deep_sky_objects[-1]

if __name__ == "__main__":
    o = OrbitalBody
    d = Convert.date2day_number()
    print Convert.day_number2date(d)
    nightb, nighte = sun.night(d)
    print 'Night:',Convert.day_number2date(nightb),
    print '-',Convert.day_number2date(nighte)
    interesting = major_solar_system + deep_sky_objects
    interesting = [pol]
    interesting = [major_solar_system[0]]
    for e in interesting:
        tra, tdec = e.equatorialt(d)
        az, alt = e.azalt(d)
        diam = e.apparent_diameter(d)

        print e.describe()
        print ' '*9,'RA/Dec %s/%s' % (Format.rad2h(tra), Format.rad2d(tdec))
        print ' '*9,'Az/Alt %s/%s' % (Format.rad2d(az), Format.rad2d(alt))
        tra2, tdec2 = Location.last.equatorial(az, alt, d)
        print ' '*9,'ra/dec %s/%s' % (Format.rad2h(tra2), Format.rad2d(tdec2))
        print ' '*9,'Mag %+6.2f / ADiam %s' % (
            e.magnitude(d), Format.rad2d(diam))
        print ' '*9,'Meridian   %s' % (
            Convert.day_number2date(e.meridian(1,d),))
        print ' '*9,'Nadir %s' % (Convert.day_number2date(e.nadir(1,d),))
        print ' '*9,'Rise %s' % (Convert.day_number2date(e.rise(-5./6,1,d),))
        print ' '*9,'Set %s' % (Convert.day_number2date(e.set(-5./6,1,d),))

        for dn,c in e.events(d, d+365, h=-50./60):
            print c, Convert.day_number2date(dn)-datetime.timedelta(hours=7)
        if not isinstance(e, OrbitalBody): continue
        elong = e.elongation(d)
        fv = e.phase_angle(d)
        dist = e.earth_distance(d)
        print ' '*9,'Elong %s / FV %s / Phase %6.4f / Dist %7.4f' % (
            Format.rad2d(elong), Format.rad2d(fv), e.phase(d), dist)
