import numpy as np
import scipy.interpolate as sci
import matplotlib.pylab as plt
import matplotlib
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyasl


class ContiInteractive:
  """
    Define the continuum estimate manually and interactively.
    
    Given wavelength and flux arrays, this tool allows to define individual
    continuum points between which either linear or spline
    interpolation is used to specify the continuum model.
    
    The main method of this class is called "searchConti". On call, it shows
    a GUI window with the spectrum. You can use the middle button of the
    mouse to *add points* defining the continuum. Once you have defined a
    sufficient number of points, the window will also show the obtained
    continuum estimate.
    
    To *remove points*, use the "remove" button. Again, clicking
    the middle button of the mouse will remove the closest point; use the "add"
    button to change the behavior again. Note that you can also use the "r"
    and "a" keys to switch between the removing and adding point modes.
    
    Additionally, there are two buttons
    entitled "cubic" and "linear", which determine whether spline or
    linear interpolation is used to calculate the continuum.
    
    Once you have obtained a reasonable estimate, you can check the result
    using the "Normalize" button. After clicking it, another window opens,
    which shows the normalized spectrum.
    
    If you are satisfied with the continuum estimate, hit the OK button to
    close the GUI window.
    
    Parameters
    ----------
    modelConti : 2d-array
        Model continuum derived from some other source. A two-dimensional
        array with the first column denoting the model wavelength and the
        second column the model flux.
    alwaysMC : boolean, optional
        If True, the model is always shown along with the data.
        Default is False.
    splineKind : string, {cubic, linear}, optional
        The type of spline to be used; can be changed via the GUI.
    useINTEP : boolean, optional
        If True, cubic interpolation will be replaced by INTEP interpolation.
  """
  
  def __init__(self, splineKind="cubic", modelConti=None, alwaysMC=False, useINTEP=False):
    self.useINTEP = useINTEP
    if self.useINTEP:
      splineKind = "intep"
    # splinekind is either "cubic" or "linear" (GUI controlled)
    self.splineKind = splineKind
    # Point list contains the entries specifying the user defined points.
    # An entry looks as follows:
    #   [xdata, ydata, True/False, 2DLine instance]
    # The third entry whether a point has already been put on the canvas
    # and the forth is the line instance, which represents that point in
    # the figure (needed to be able to remove it).
    self.pointList = []
    # The current mode ("add" points or "remove" points)
    self.mode = "add"
    # Continuum model
    self.cmodel = None
    # Connection of middle mouse button
    self.middleButtonCon = None
    # The model continuum
    self.modelConti = modelConti
    # Always show model continuum?
    self.alwaysMC = alwaysMC
    # Saves instance of the line used to plot comparison model 
    self.modelMC = None
    if self.alwaysMC and (modelConti is None):
      raise(PE.PyAParameterConflict("If alwaysMC is set True, modelConti must be provided, too."))
    
    # Inactive Button Color
    self.ibc = [0.7,0.7,0.7]
    # Active Button Color
    self.abc = [0.1,0.1,0.99]
  
  def __comparePoints(self, a, b):
    """
      Compare two points from the point list according to xadata position.
      Used to sort the point list.
    """
    if a[0] < b[0]: return -1
    if a[0] == b[0]: return 0
    return 1
    
  def getModel(self):
    """
      Get the currently defined continuum model.
      
      Either linear or spline interpolation is used to connect the
      individual points.
      
      Returns
      -------
      Continuum : array
          The continuum model (note that it
          contains `Nan' values where no continuum is defined).
    """
    if self.splineKind == "cubic":
      if len(self.pointList) < 4: return None
    if (self.splineKind == "linear") or (self.splineKind == "intep"):
      if len(self.pointList) < 2: return None
    self.__sortPointList()
    x = []; y = []
    for p in self.pointList:
      x.append(p[0]); y.append(p[1])
    if self.splineKind == "intep":
      return pyasl.intep(np.array(x), np.array(y), self.w, boundsError=False, fillValue=np.NaN)
    ip = sci.interp1d(x, y, kind=self.splineKind, bounds_error=False)
    return ip(self.w)
    
  def __sortPointList(self):
    self.pointList.sort(self.__comparePoints)
  
  def __plot(self):
    """
      Put points and continuum model into the figure. To do so, the old
      continuum model has to be removed.
    """
    self.__sortPointList()
    # Plot the points (if not already in the figure).
    for i in xrange(len(self.pointList)):
      if not self.pointList[i][2]:
        # Show only those points, which are not too far off the
        # visible axis
        range = self.w.max() - self.w.min()
        if (self.pointList[i][0] < self.w.min()-0.1*range) or \
           (self.pointList[i][0] > self.w.max()+0.1*range): continue 
        self.ax.plot([self.pointList[i][0]], [self.pointList[i][1]], 'rp')[0]
        self.pointList[i][2] = True
        self.pointList[i][3] = self.ax.lines[-1]
    
    # Plot the model onto the canvas and remove the old one
    # if necessary.
    model = self.getModel()
    if self.cmodel is None:
      if model is not None:
        self.ax.plot(self.w, self.getModel(), 'r--')
        self.cmodel = self.ax.lines[-1]
    else:
      for i in xrange(len(self.ax.lines)):
        if self.ax.lines[i] is self.cmodel:
          self.ax.lines.pop(i)
          self.cmodel = None
          break
      if model is not None:
        self.ax.plot(self.w, self.getModel(), 'r--')
        self.cmodel = self.ax.lines[-1]
    # Plot the comparison model
    if self.alwaysMC:
      if not self.modelMC is None:
        for i in xrange(len(self.ax.lines)):
          if self.ax.lines[i] is self.modelMC:
            self.ax.lines.pop(i)
            self.modelMC = None
          break
      self.ax.plot(self.modelConti[::,0], self.modelConti[::,1], 'g--')
      self.modelMC = self.ax.lines[-1]
  
  def __keyEvent(self, event):
    """
      The key ``a'' can be used to change mode to add points, while ``r''
      will start the remove point mode.
    """
    if event.key == "a":
      self.mode = "add"
    elif event.key == "r":
      self.mode = "remove"
    self.__setConnections()
  
  def __addPoint(self, event):
    """
      Add a point to the list and replot.
    """
    # Accept only middle button.
    if event.button != 2: return
    if (event.xdata is None) or (event.ydata is None):
      return
    self.pointList.append([event.xdata, event.ydata, False, None])
    self.__plot()
    self.fig.canvas.draw()
  
  def __removePoint(self, event):
    """
      Remove a point from the list and replot.
    """
    if len(self.pointList) == 0: return
    # Check whether click occurred on plot axes.
    if event.inaxes is not self.ax: return
    # Accept only middle button
    if event.button != 2: return
    # Find the index of the point closest to the click.
    dist = []
    for p in self.pointList:
      dist.append(np.sqrt((p[0]-event.xdata)**2 + (p[1]-event.ydata)**2))
    dist = np.array(dist)
    indmin = np.argmin(dist)
    # Remove that point and redraw.
    for i in xrange(len(self.ax.lines)):
      if self.ax.lines[i] is self.pointList[indmin][3]:
        self.ax.lines.pop(i)
        self.pointList.pop(indmin)
        break
    self.__plot()
    self.fig.canvas.draw()
    
  
  def __setConnections(self):
    """
      Depending on the mode (add or remove) change meaning of the middle
      mouse button.
    """
    if self.mode == "add":
      if self.middleButtonCon is not None:
        self.fig.canvas.mpl_disconnect(self.middleButtonCon)
      self.middleButtonCon = self.fig.canvas.mpl_connect("button_press_event", self.__addPoint)
    if self.mode == "remove":
      if self.middleButtonCon is not None:
        self.fig.canvas.mpl_disconnect(self.middleButtonCon)
      self.middleButtonCon = self.fig.canvas.mpl_connect("button_press_event", self.__removePoint)
  
  def __cubicClicked(self, event):
    """
      The cubic button click event handler. Change splineKind and button
      color.
    """
    if not self.useINTEP:
      self.splineKind = "cubic"
    else:
      self.splineKind = "intep"
    # Change button color
    self.buttonCubic.color = self.abc
    self.buttonLinear.color = self.ibc
    # Plot (note that _motio calls are needed to change
    # the color immediately on click).
    self.__plot()
    self.buttonCubic._motion(event)
    self.buttonLinear._motion(event)
    self.fig.canvas.draw()

  def __linearClicked(self, event):
    self.splineKind = "linear"
    # Change button color
    self.buttonCubic.color = self.ibc
    self.buttonLinear.color = self.abc
    # Plot (note that _motio calls are needed to change
    # the color immediately on click).
    self.__plot()
    self.buttonCubic._motion(event)
    self.buttonLinear._motion(event)
    self.fig.canvas.draw()

  def __addClicked(self, event):
    """
      Change mode to adding points.
    """
    self.mode = "add"
    # Plot (note that _motio calls are needed to change
    # the color immediately on click).
    self.buttonAdd.color = self.abc
    self.buttonRem.color = self.ibc
    self.buttonAdd._motion(event)
    self.buttonRem._motion(event)
    self.__setConnections()

  def __remClicked(self, event):
    """
      Change mode to remove points.
    """
    self.mode = "remove"
    # Plot (note that _motio calls are needed to change
    # the color immediately on click).
    self.buttonAdd.color = self.ibc
    self.buttonRem.color = self.abc
    self.buttonAdd._motion(event)
    self.buttonRem._motion(event)
    self.__setConnections()

  def __normClicked(self, event):
    """
      Normalization button was clicked. Open a second figure and show
      the normalized spectrum to check the result.
    """
    fign = plt.figure(facecolor="white")
    ax = fign.add_subplot(111, title="Check normalization")
    model = self.getModel()
    nflux = self.f / model
    nnone = np.where(np.isnan(model) == False)[0]
    ax.plot(self.w[nnone], nflux[nnone],  'b-')
    ax.plot([self.w[nnone].min(), self.w[nnone].max()], [1,1], 'r--')
    if self.modelConti is not None:
      indi = np.where(np.logical_and( \
             self.modelConti[::,0] > self.w[nnone].min(),\
             self.modelConti[::,0] < self.w[nnone].max()))[0]
      ax.plot(self.modelConti[indi,0], self.modelConti[indi,1], 'g--')
    fign.show()

  def __okClicked(self, event):
    plt.close()

  def __defineButtons(self):
    """
      Put buttons on the canvas and assign meaning.
    """
    # Linear and cubic
    axButtonCubic = self.fig.add_axes([0.1, 0.0375, 0.1, 0.05])
    axButtonLinear = self.fig.add_axes([0.22, 0.0375, 0.1, 0.05])
    
    # CubicButton text
    cubt = "cubic"
    if self.useINTEP:
      # Replace of INTEP
      cubt = "INTEP"
    
    if (self.splineKind == "cubic") or (self.splineKind == "intep"):
      self.buttonCubic = matplotlib.widgets.Button(axButtonCubic, cubt, color=self.abc)
      self.buttonLinear = matplotlib.widgets.Button(axButtonLinear, "linear", color=self.ibc)
    else:
      self.buttonCubic = matplotlib.widgets.Button(axButtonCubic, cubt, color=self.ibc)
      self.buttonLinear = matplotlib.widgets.Button(axButtonLinear, "linear", color=self.abc)
    self.buttonCubic.on_clicked(self.__cubicClicked)
    self.buttonLinear.on_clicked(self.__linearClicked)
    
    # Add and remove
    axButtonAdd = self.fig.add_axes([0.36, 0.0375, 0.1, 0.05])
    axButtonRem = self.fig.add_axes([0.48, 0.0375, 0.1, 0.05])
    self.buttonAdd = matplotlib.widgets.Button(axButtonAdd, "add", color=self.abc)
    self.buttonRem = matplotlib.widgets.Button(axButtonRem, "rem", color=self.ibc)  
    self.buttonAdd.on_clicked(self.__addClicked)
    self.buttonRem.on_clicked(self.__remClicked)
    
    # Normalize
    axbuttonNorm = self.fig.add_axes([0.65, 0.0375, 0.15, 0.05])
    self.buttonNorm = matplotlib.widgets.Button(axbuttonNorm, "Normalize")
    self.buttonNorm.on_clicked(self.__normClicked)
    
    # OK
    axButtonOK = self.fig.add_axes([0.83, 0.0375, 0.07, 0.05])
    self.buttonOK = matplotlib.widgets.Button(axButtonOK, "OK")
    self.buttonOK.on_clicked(self.__okClicked)

  def __cleanPointListAfterFigureClose(self):
    """
      If the figure is closed, the plot flag and the line reference entries \
      in the pointList loose their meaning and have to be reset. In this way, \
      the points are again available, when the plot is created (searchConti \
      called) the next time.
    """
    for i in xrange(len(self.pointList)):
      self.pointList[i][2] = False
      self.pointList[i][3] = None

  def calc(self, w, f, points=[], region=1., interpol = 'cubic'):
      """
        Non-interactive method to determine the continuum by estimating the 
        continuum level at the given points.
        
        Parameters
        ----------
          w : array 
               The wavelength array
          f : array
               The flux array
          points : array or list
               The wavelength points where the continuum will be estimated
          region : float or array
               Wavelength range around the given points for continuum estimate
          interpol : string
               Defines the interpolation scheme, must be either 'linear' or 'cubic'                        
      """
      self.w = w.copy()
      self.f = f.copy()
      self.pointList = []
      self.splineKind = interpol
      if type(region) in [type(0.1),type(1)]:
          region = np.ones(len(points))*region
      for i in range(len(points)):
          gindex = np.where(abs(self.w-points[i]) < region[i]/2.)[0]
          y = np.median(self.f[gindex])
          self.pointList.append([points[i], y, False, None])   
      m =  self.getModel()
      return m
  
  def searchConti(self, w, f, figTitle="", compare=None):
    """
      Find the continuum estimate manually and interactively.
      
      Parameters
      ----------
      w : array
          The wavelength array
      f : array
          The flux array
      figTitle : string, optional
          Defines the title of the figure.
      compare : array
          An array with two columns: wvl and flux, which will
          be plotted along the data for comparison in the
          window used to define the continuum.
      
      Returns
      -------
      Continuum : array
          The continuum model evaluated at the input wavelength
          points.
    """
    
    # Create a figure
    self.fig = plt.figure(facecolor="white")
    # Handle key pressed events
    self.fig.canvas.mpl_connect("key_press_event", self.__keyEvent)
    self.ax = self.fig.add_subplot(111)
    # Leave a little room for buttons
    self.fig.subplots_adjust(bottom=0.15)

    # Put the buttons into the figure
    self.__defineButtons()
    
    self.ax.set_title(figTitle)
    self.w = w.copy()
    self.f = f.copy()
    
    self.ax.plot(self.w, self.f, 'b--')
    if compare is not None:
      self.ax.plot(compare[::,0], compare[::,1], 'g--')
    
    self.__plot()
    self.__setConnections()

    plt.show()
    self.__cleanPointListAfterFigureClose()
    return self.getModel()
  
  def __call__(self, w, f, **kwargs):
    return self.searchConti(w, f, "Manual continuum approximation", **kwargs)
