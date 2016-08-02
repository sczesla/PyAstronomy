# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy
import os
from . import astroTimeLegacy as at
from . import aitoffLegacy
import unittest
from .baryvel import baryvel, helcorr
from . import eq2hor
from . import sunpos
from . import moonpos
from . import moonphase
from . import posAngle
from .airtovac import airtovac2, vactoair2
import six.moves as smo

#  How to run...
#  p3 -m PyAstronomy.pyasl.asl.test
#

class AstroTimeLegacyTest(unittest.TestCase):

  def setUp(self, tddir="testPro", p=1e-6):
    """
      Parameter:
        tddir - Directory where test data can be found.
    """
    self.tddir = tddir
    self.p = p

  def getData(self, fn):
    if not os.path.isfile(self.tddir+"/"+fn):
      print("Could not find test data file: ", fn)
      print("  Use test.pro (create_test_data) to generate test data.")
      return [False, None]
    return [True, numpy.loadtxt(self.tddir+"/"+fn)]

  def test_bprecess(self):
    fn = "bprecess.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    
    err = 0
    for i in smo.range(len(dat[::,0])):
      ra, de, r5, d5, pa1, pa, m1_1, m2_1, m1_2, m2_2, rv1, rv = dat[i,::]
      mu = [m1_1, m2_1]
      if m1_1 == 0.0 and m2_1 == 0.0:
        mu = None
      if pa1 == 0.0:
        pa = 0.0
      if rv1 == 0.0:
        rv = 0.0
      
      r = at.bprecess(ra, de, parallax=pa, mu_radec=mu, rad_vel=rv)
      # [ra_1950, dec_1950, MU_RADEC, PARALLAX, RAD_VEL]
      
      if (r5-r[0])/r5 > self.p:
        print("Problem r5: ", (r5-r[0])/r5)
        err += 1
        print(dat[i,::])
        print(r)
        print("")
        continue
      if (d5-r[1])/d5 > self.p:
        print("Problem d5: ", (d5-r[1])/d5)
        err += 1
        print(dat[i,::])
        print(r)
        print("")
        continue
      if r[2] is not None:
        if (r[2][0]-m1_2)/m1_2 > self.p:
          print("Problem m1: ", (r[2][0]-m1_2)/m1_2)
          err += 1
          print(dat[i,::])
          print(r)
          print("")
          continue
        if (r[2][1]-m2_2)/m2_2 > self.p:
          print("Problem m2: ", (r[2][1]-m2_2)/m2_2)
          err += 1
          print(dat[i,::])
          print(r)
          print("")
          continue
      if r[3] != 0.0:
        if (r[3]-pa)/pa > self.p:
          print("Problem parallax: ", (r[3]-pa)/pa)
          err += 1
          print(dat[i,::])
          print(r)
          print("")
          continue
      if rv != 0.0:
        if (r[4]-rv)/rv > self.p:
          print("Problem rv: ", (r[4]-rv)/rv)
          err += 1
          print(dat[i,::])
          print(r)
          print("")
          continue
    
    self.assertEqual(err, 0)

  def test_premat(self):
    fn = "premat.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      r = numpy.reshape(dat[i,2:11], (3,3))
      if dat[i,11] == 0:
        m = at.premat(dat[i,0], dat[i,1])
      else:
        m = at.premat(dat[i,0], dat[i,1], FK4=True)
      if (r-m).max() > self.p:
        print("test_premat - Error detected")
        err += 1
    self.assertEqual(err, 0)

  def test_precess(self):
    fn = "precess.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      fk4 = True
      rad = True
      if dat[i,6] == 0.0:
        fk4 = False
      if dat[i,7] == 0.0:
        rad = False
      r = at.precess(dat[i,2], dat[i,3], dat[i,0], dat[i,1], FK4=fk4, radian=rad)
      if (r[0]-dat[i,4])/r[0] > self.p:
        err += 1
      if (r[1]-dat[i,5])/r[0] > self.p:
        err += 1
    self.assertEqual(err, 0)


  def test_precess_xyz(self):
    fn = "precess_xyz.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      x,y,z = at.precess_xyz(dat[i,2], dat[i,3], dat[i,4], dat[i,0], dat[i,1])
      if abs((x-dat[i,5])/x) > self.p:
        err += 1
        print("Probelm x: ", (x-dat[i,5])/x)
        print(x,y,z)
        print(dat[i,5:8])
        print(numpy.array([x,y,z]) - dat[i,5:8])
        print("")
        continue
      if abs((y-dat[i,6])/y) > self.p:
        err += 1
        print("Probelm y: ", (y-dat[i,6])/y)
        print(x,y,z)
        print(dat[i,5:8])
        print(numpy.array([x,y,z]) - dat[i,5:8])
        print("")
        continue
      if abs((z-dat[i,7])/z) > self.p:
        err += 1
        print("Probelm z: ", (z-dat[i,7])/z)
        print(x,y,z)
        print(dat[i,5:8])
        print(numpy.array([x,y,z]) - dat[i,5:8])
        print("")
        continue
    self.assertEqual(err, 0)

  def test_xyz(self):
    fn = "xyz.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      r = numpy.array(at.xyz(dat[i,0], velocity=True, equinox=dat[i,7]))
      m = (numpy.abs(r-dat[i,1:7])/r).max()
      if m > self.p:
        err += 1
        print("Problem: ", m)
        print(r)
        print(m)
        print("")
        continue
    self.assertEqual(err, 0)

  def test_helio_jd(self):
    fn = "helio_jd.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      hjd = at.helio_jd(dat[i,1], dat[i,2], dat[i,3])
      if abs((hjd-dat[i,0])/hjd) > self.p:
        print("Problem HJD: ", abs((hjd-dat[i,0])/hjd))
        print(hjd, dat[i,0])
        print("")
        err += 1
        continue
    self.assertEqual(err, 0)

  def test_daycnv(self):
    fn = "daycnv.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      r = numpy.array(at.daycnv(dat[i,0]))
      m = numpy.abs((r - dat[i,1::])/r).max()
      if m > self.p:
        print("Problem: ")
        print(r)
        print(dat[i,1::])
        err += 1
        continue
    self.assertEqual(err, 0)
  
  def test_aitoff(self):
    fn = "aitoff.test"
    good, dat = self.getData(fn)
    self.assertEqual(good, True)
    err = 0
    for i in smo.range(len(dat[::,0])):
      x, y = aitoffLegacy.aitoff(dat[i,0], dat[i,1])
      # The second condition (abs(x)>1e-2)) because IDL uses single precision
      if (abs((x-dat[i,2])/x) > 1e-4 and (abs(x)>1e-2)) or (abs((y-dat[i,3])/y) > 1e-4 and (abs(y) > 1e-2)):
        print("Problem: ")
        print("  l, b, x(IDL), y(IDL), x(Python), y(Python): ", dat[i,0], dat[i,1], x, y, dat[i,2], dat[i,3])
        err += 1
    self.assertEqual(err, 0)

  def test_inverseAitoff(self):
    err = 0
    for i in smo.range(1000):
      l = numpy.random.random() * 360.0 - 180.0
      b = numpy.random.random() * 180.0 - 90.0
      x, y = aitoffLegacy.aitoff(l, b)
      l2, b2 = aitoffLegacy.inverseAitoff(x, y)
      if abs(l-l2) > 1e-6 or abs(b-b2) > 1e-6:
        print("Problem: ")
        print("  ", l, b, x, y, l2, b2)
        err += 1
    self.assertEqual(err, 0)
  
  def test_baryvel(self):
    dat = numpy.loadtxt("testPro/baryvel.test")
    err = 0
    for i in range(10000):
      jd = 2.4e6 + float(i)*10.37681
      a, b =  baryvel(jd, 0)
      delta = numpy.abs(dat[i,::] - numpy.concatenate((a,b)))
      if max(delta) > 1e-6:
        err += 1
        print("Baryvel error:", dat[i,::], a, b)
        print("max: ", max(delta))
    self.assertEqual(err, 0)

    i = 0
    dat = numpy.loadtxt("testPro/baryvel2.test")
    for j in range(101):
      jd = 2.4e6 + float(j)*1000.37681
      for k in range(11):
        deq = jd + float(k-5) * 10000.0
        a, b =  baryvel(jd, deq)
        delta = numpy.abs(dat[i,0:6] - numpy.concatenate((a,b)))
        
        if abs(jd - dat[i,6]) > 1e-6:
          print("JD problem: ", jd, dat[i,6], i, k)
        if abs(deq - dat[i,7]) > 1e-6:
          print("deq problem: ", deq, dat[i,7], i, k)
        
        if max(delta) > 1e-6:
          err += 1
          print("Baryvel error:", dat[i,::], a, b)
          print("max: ", max(delta))
          print("i, k: ", i, k)
        i += 1


class IDLTests(unittest.TestCase):

  def setUp(self, tddir="testPro", p=1e-6):
    """
      Parameter:
        tddir - Directory where test data can be found.
    """
    self.tddir = tddir
    self.p = p

  def getData(self, fn):
    if not os.path.isfile(self.tddir+"/"+fn):
      print("Could not find test data file: ", fn)
      print("  Use test.pro (create_test_data) to generate test data.")
      return [False, None]
    return [True, numpy.loadtxt(self.tddir+"/"+fn)]
    
  def test_sunpos(self):
    """
      Testing sunpos
    """
    dat = self.getData("sunpos.test")[1]
    for i in smo.range(len(dat[::,0])):
      jd, ra, dec, longmed, oblt = sunpos.sunpos(dat[i,0], full_output=True)
      self.assertAlmostEqual(ra, dat[i,1], delta=self.p)
      self.assertAlmostEqual(dec, dat[i,2], delta=self.p)
      self.assertAlmostEqual(longmed, dat[i,3], delta=self.p)
      self.assertAlmostEqual(oblt, dat[i,4], delta=self.p)

  def test_nutate(self):
    dat = self.getData("nutate.test")[1]
    for i in smo.range(len(dat[::,0])):
      l, o = eq2hor.nutate(dat[i,0])
      self.assertAlmostEqual(l*3600.0, dat[i,1], delta=self.p)
      self.assertAlmostEqual(o*3600.0, dat[i,2], delta=self.p)   
    
  def test_co_nutate(self):
    dat = self.getData("co_nutate.test")[1]
    for i in smo.range(len(dat[::,0])):
      dra, ddec, o, dl, do = eq2hor.co_nutate(dat[i,0], dat[i,1], dat[i,2], full_output=True)
      self.assertAlmostEqual(dra*3600./dat[i,3], 1.0, delta=1e-5)
      self.assertAlmostEqual(ddec*3600./dat[i,4], 1.0, delta=1e-5)
      self.assertAlmostEqual(o/180.*numpy.pi/dat[i,5], 1.0, delta=self.p)
      self.assertAlmostEqual(dl*3600./dat[i,6], 1.0, delta=self.p)
      self.assertAlmostEqual(do*3600./dat[i,7], 1.0, delta=self.p)

  def test_co_aberration(self):
    dat = self.getData("co_aberration.test")[1]
    for i in smo.range(len(dat[::,0])):
      dra, ddec = eq2hor.co_aberration(dat[i,0], dat[i,1], dat[i,2])
      self.assertAlmostEqual(dra*3600./dat[i,3], 1.0, delta=self.p)
      self.assertAlmostEqual(ddec*3600./dat[i,4], 1.0, delta=self.p)

  def test_co_refract_forward(self):
    dat = self.getData("co_refract_forward.test")[1]
    for i in smo.range(len(dat[::,0])):
      r = eq2hor.co_refract_forward(dat[i,0], pressure=dat[i,1], temperature=dat[i,2]-273.15)
      self.assertAlmostEqual(r,dat[i,3], delta=self.p)

  def test_co_refract(self):
    dat = self.getData("co_refract.test")[1]
    for i in smo.range(len(dat[::,0])):
      aout, p, t = eq2hor.co_refract(dat[i,0], observer_alt=dat[i,1], pressure=dat[i,2], \
                               temperature=dat[i,3]-273.15, convert_to_observed=True)
      self.assertAlmostEqual(aout/dat[i,4], 1.0, delta=self.p)

  def test_hadec2altaz(self):
    dat = self.getData("hadec2altaz.test")[1]
    for i in smo.range(len(dat[::,0])):
      alt, az = eq2hor.hadec2altaz(dat[i,0], dat[i,1], dat[i,2])
      self.assertAlmostEqual(alt/dat[i,3], 1.0, delta=self.p)
      self.assertAlmostEqual(az/dat[i,4], 1.0, delta=self.p)

  def test_eq2hor(self):
    """
      Testing eq2hor
    """
    dat = self.getData("eq2hor.test")[1]
    for i in smo.range(len(dat[::,0])):
      alt, az, ha = eq2hor.eq2hor(dat[i,2], dat[i,0], dat[i,1], lon=dat[i,3], lat=dat[i,4], alt=dat[i,5])
      self.assertAlmostEqual(alt/dat[i,6], 1.0, delta=self.p*100., msg="Altitude does not fit " + str(dat[i,::]))
      self.assertAlmostEqual(az/dat[i,7], 1.0, delta=self.p*100., msg="Azimuth does not fit " + str(dat[i,::]))
  
  def test_moonpos(self):
    """
      Testing moonpos routine
    """
    dat = self.getData("moonpos.test")[1]
    for i in smo.range(len(dat[::,0])):
      ra, dec, dist, glon, glat = moonpos.moonpos(dat[i,0], False)
      self.assertAlmostEqual(ra/dat[i,1], 1.0, delta=self.p)
      self.assertAlmostEqual(dec/dat[i,2], 1.0, delta=self.p)
      self.assertAlmostEqual(dist/dat[i,3], 1.0, delta=self.p)
      self.assertAlmostEqual(glon/dat[i,4], 1.0, delta=self.p)
      self.assertAlmostEqual(glat/dat[i,5], 1.0, delta=self.p)

  def test_mphase(self):
    """
      Testing moonphase (mphase) routine
    """
    dat = self.getData("mphase.test")[1]
    for i in smo.range(len(dat[::,0])):
      f = moonphase.moonphase(dat[i,0])
      self.assertAlmostEqual(f, dat[i,1], delta=1e-5, msg="Lunar phase does not match." )
      self.assertAlmostEqual(f/dat[i,1], 1.0, delta=1e-5, msg="Lunar phase (relative) does not match")

  def test_posangle(self):
    """
      Testing positionAngle (posAng) routine
    """
    dat = self.getData("posangle.test")[1]
    for i in smo.range(len(dat[::,0])):
      # Using *15.0 because decimal hours have been specified for IDL...
      f = posAngle.positionAngle(dat[i,0]*15.0, dat[i,1], dat[i,2]*15.0, dat[i,3], positive=False)
      self.assertAlmostEqual(dat[i,4], f, delta=1e-9, msg="Position angle does not match.")

  def test_helcorr(self):
    """
      Testing helcorr function for barycentric correction
    """
    dat = self.getData("helcorr.test")[1]
#    import matplotlib.pylab as plt
#    diffs = []
    for i in smo.range(len(dat[::,0])):
      # Using *15.0 because decimal hours have been specified for IDL...
      # The minus accounts for the east-west change
      corr, hjd = helcorr(-dat[i,0], dat[i,1], dat[i,2], dat[i,3]*15.0, dat[i,4], dat[i,5])
#      print corr, dat[i,6], corr-dat[i,6]
#      diffs.append(corr-dat[i,6])
      self.assertAlmostEqual(dat[i,6], corr, delta=1e-5, msg="Barycentric correction does not match.")
      self.assertAlmostEqual(dat[i,7]+2.4e6, hjd, delta=1e-5, msg="HJD does not match.")
#    plt.hist(diffs)
#    plt.show()

  def test_airtovac(self):
    """
      Testing air to vac conversion
    """
    dat = self.getData("airvac2.test")[1]
    wvl = dat[::,0]
    atv = dat[::,1]
    vta = dat[::,2]
    atv2 = airtovac2(wvl)
    vta2 = vactoair2(wvl)
    print("test_airtovac")
    print("  max diff in airtovac", numpy.max(numpy.abs(atv-atv2)))
    print("  max diff in vactoair", numpy.max(numpy.abs(vta-vta2)))
    self.assertAlmostEqual(numpy.max(numpy.abs(atv-atv2)), 0.0, delta=1e-6, msg="Air to Vac conversion inconsistent.")
    self.assertAlmostEqual(numpy.max(numpy.abs(vta-vta2)), 0.0, delta=1e-6, msg="Vac to air conversion inconsistent.")

if __name__ == "__main__":
  suite = unittest.TestLoader().loadTestsFromTestCase(AstroTimeLegacyTest)
  unittest.TextTestRunner(verbosity=2).run(suite)
  
  suite2 = unittest.TestLoader().loadTestsFromTestCase(IDLTests)
  unittest.TextTestRunner(verbosity=2).run(suite2)
