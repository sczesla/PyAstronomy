import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf
import numpy as np

import six.moves.tkinter as tk
import six.moves as smo
import six


class Point():
  
  def __init__(self, event=None):
    self.xdata = 0.0
    self.ydata = 0.0
    self.mplLine = None
    self.lbIdent = None
    if event is not None:
      self.xdata = event.xdata
      self.ydata = event.ydata
  
  def asStr(self):
    """
      String representation of the point.
      
      Returns
      -------
      Point : string
          String representation of the point coordinates.
    """
    return "%g, %g" % (self.xdata, self.ydata)


class IAGVFit:
  """
    GUI for an interactive of spectral lines.
    
    Parameters
    ----------
    x, y : array
        The data (e.g., spectrum).
    yerr : array, optional
        The error on the data.
    mode : string, {voigt,gauss}, optional
        The line profile to be used. Either a
        Gaussian or a Voigt profile can be used.
    config : dictionary, optional
        A dictionary containing any combination of
        the keys, which can be used to change the default
        behavior of the class.

        +------------------------+------------+---------------------------------+
        | Key                    | Default    | Description                     |
        +========================+============+=================================+
        | modelLS                | r--        | Style used for the combined     |
        |                        |            | model.                          |
        +------------------------+------------+---------------------------------+
        | dataLS                 | bp         | Style used for the data of no   |
        |                        |            | error is given.                 |
        +------------------------+------------+---------------------------------+
        | dataLSerr              | b+         | Style used when error is given. |
        +------------------------+------------+---------------------------------+
        | activeCompLS           | k-         | Style used to plot active       |
        |                        |            | component.                      |
        +------------------------+------------+---------------------------------+
        | inactiveCompLS         | y--        | Style used to plot inactive     |
        |                        |            | component(s).                   |
        +------------------------+------------+---------------------------------+
        | plotIndComps           | True       | If True, individual components  |
        |                        |            | are plotted.                    |
        +------------------------+------------+---------------------------------+
        | rangeLS                | m--        | Style used to plot lines        |
        |                        |            | indicating the fitting range.   |
        +------------------------+------------+---------------------------------+
        | defaultFreeVoigt       | A, ad, al, | Default free fitting parameters |
        |                        | mu         | for Voigt profile.              |
        +------------------------+------------+---------------------------------+
        | default_alad           | 0.1        | Default ratio between width of  |
        |                        |            | Lorentz and Gauss contribution. |
        +------------------------+------------+---------------------------------+
        | defaultFreeGauss       | A, ad, sig | Default set of free fitting     |
        |                        |            | parameters for Gaussian.        |
        +------------------------+------------+---------------------------------+
        | defaultWheelAdd        | 0.01       | Default increment for additive  |
        |                        |            | mouse wheel manipulation.       |
        +------------------------+------------+---------------------------------+
        | defaultWheelMult       | 2          | Default factor for              |
        |                        |            | multiplicative mouse wheel      |
        |                        |            | manipulation.                   |
        +------------------------+------------+---------------------------------+
                    
    Attributes
    ----------
    f : mpl Figure
        A Figure instance from matplotlib.
    a : mpl axis
        An Axis instance from matplotlib.
  """
  
  def __init__(self, x, y, yerr=None, mode="gauss", config=None):
    self.windowTitle = "PyA interactive GV"
    self.f = Figure()
    self.a = self.f.add_subplot(111)
        
    self._config = {"modelLS":'r--', \
                    "dataLS":'bp', \
                    "dataLSerr":'b+', \
                    "activeCompLS":'k-', \
                    "inactiveCompLS":'y--', \
                    "plotIndComps":True, \
                    "rangeLS":'m--', \
                    "defaultFreeVoigt":["A", "ad", "al", "mu"], \
                    "default_alad":0.1, \
                    "defaultFreeGauss":["A", "sig", "mu"], \
                    "defaultWheelAdd":0.01, \
                    "defaultWheelMult": 2}
    
    # Modify configuration according to input
    if config is not None:
      for k in six.iterkeys(self._config):
        if k in config:
          # Value has been modified
          self._config[k] = config[k]
    
    # Tk root
    self.root = tk.Tk()
    
    # Create data plot
    self._x = x
    self._y = y
    self._yerr = yerr
    if self._yerr is None:
      self.a.plot(x, y, self._config["dataLS"])
    else:
      self.a.errorbar(x, y, yerr=self._yerr, fmt=self._config["dataLSerr"])
    
    # The input mode determines how clicks on the mouse wheel are
    # interpreted.
    self._inputMode = "newComp"
     
    # Component count/assignment starts at 1
    self._activeComponent = 1
    self._compoCount = 0
    
    # Create model and parameter lists
    self._mode = mode
    if self._mode == "voigt":
      self._model = fuf.MultiVoigt1d(1)
      self._compPars = ["mu", "ad", "al", "A"]
    elif self._mode == "gauss":
      self._model = fuf.MultiGauss1d(1)
      self._compPars = ["mu", "sig", "A"]
    else:
      raise(PE.PyAValError("No such mode: " + str(self._mode), \
                           solution="Choose either 'voigt' or 'gauss'."))
    self._model["off"] = 1.0

    # List used to cache points used to add component
    # as long as component is not completely specified
    self._pointCache = []

    # What to consider in the fit
    self._rangeIndices = np.arange(self._x.size, dtype=np.int)
    
    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.pack(fill=tk.BOTH, side=tk.LEFT, expand=True)
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.pointFrame = tk.Frame(self.root)
    self.pointFrame.pack(side=tk.LEFT, fill=tk.BOTH)
    
    # Manage active component
    self.activeCompoFrame = tk.Frame(self.pointFrame)
    self.activeCompoFrame.pack(side=tk.TOP)
    self.currentComponentLabel = tk.Label(self.activeCompoFrame, text="Active component")
    self.currentComponentLabel.pack(side=tk.TOP)
    # De- and increase no. of active component
    self.downACB = tk.Button(self.activeCompoFrame, text="-", command=self._downAC)
    self._acstring = tk.StringVar()
    self._acstring.set(str(self._activeComponent)) 
    self.ACLabel = tk.Label(self.activeCompoFrame, textvariable=self._acstring)
    self.upACB = tk.Button(self.activeCompoFrame, text="+", command=self._upAC)
    self.downACB.pack(side=tk.LEFT)
    self.ACLabel.pack(side=tk.LEFT)
    self.upACB.pack(side=tk.RIGHT)
    
    # Show and work with parameters
    self.parameterFrame = tk.Frame(self.pointFrame, relief=tk.RAISED, borderwidth=1)
    self.parameterFrame.pack(side=tk.TOP)
    # Headline
    tk.Label(self.parameterFrame, text="Parameters (active comp.)").pack(side=tk.TOP)
    # List of parameters to show
    pars = ["off"]
    pars.extend(self._compPars)
    self._parSelect = tk.IntVar()
    # Frames for individual parameters
    self._parFrames = []
    self._parValues = []
    self._parFreeVars = []
    
    for i, p in enumerate(pars):
      self._parFrames.append(tk.Frame(self.parameterFrame))
      self._parFrames[-1].pack(side=tk.TOP)
      b = tk.Radiobutton(self._parFrames[-1], text="", variable=self._parSelect, value=i)
      b.pack(side=tk.LEFT, anchor=tk.W)
      l = tk.Label(self._parFrames[-1], text=("%3s" % p), width=3)
      l.pack(side=tk.LEFT, anchor=tk.W)
      self._parValues.append(tk.StringVar())
      self._parValues[-1].set("0.0")
      l = tk.Label(self._parFrames[-1], textvariable=self._parValues[-1], width=10)
      l.pack(side=tk.LEFT)
      
      self._parFreeVars.append(tk.IntVar())
      c = tk.Checkbutton(self._parFrames[-1], text="free", variable=self._parFreeVars[-1], \
                         command=self._freeCheckboxClicked)
      c.pack(side=tk.LEFT)
    
    # Wheel modification
    self._wheelModFrame = tk.Frame(self.pointFrame, relief=tk.RAISED, borderwidth=1)
    self._wheelModFrame.pack(side=tk.TOP)
    # Headline
    tk.Label(self._wheelModFrame, text="Mouse wheel").pack(side=tk.TOP)
    self._multAddSelect = tk.IntVar()
    self._multAddSelect.set(1)
    multiFrame = tk.Frame(self._wheelModFrame)
    addFrame = tk.Frame(self._wheelModFrame)
    multiFrame.pack(side=tk.TOP)
    addFrame.pack(side=tk.TOP)
    tk.Radiobutton(multiFrame, text="multiply [%]", variable=self._multAddSelect, value=1, width=9, anchor=tk.W).pack(side=tk.LEFT)
    tk.Radiobutton(addFrame, text="add", variable=self._multAddSelect, value=2, width=9, anchor=tk.W).pack(side=tk.LEFT)
    self._multNumber = tk.StringVar()
    self._addNumber = tk.StringVar()
    self._multNumber.set(str(self._config["defaultWheelMult"]))
    self._addNumber.set(str(self._config["defaultWheelAdd"]))
    self._multNumberEntry = tk.Entry(multiFrame, textvariable=self._multNumber, width=10)
    self._multNumberEntry.pack(side=tk.RIGHT)
    self._addNumberEntry = tk.Entry(addFrame, textvariable=self._addNumber, width=10)
    self._addNumberEntry.pack(side=tk.RIGHT)
    
    # Common width
    if self._mode == "gauss":
      commonSigFrame = tk.Frame(self.pointFrame)
      commonSigFrame.pack(side=tk.TOP)
      self._commonSig = tk.IntVar()
      c = tk.Checkbutton(commonSigFrame, text="Common width (sig)", variable=self._commonSig, \
                         command=self._commonSigClicked)
      c.pack(side=tk.LEFT)
    
    # Fit button
    self._fitButton = tk.Button(self.pointFrame, text="Fit", command=self._fitModel)
    self._fitButton.pack(side=tk.TOP, fill=tk.X)
    
    # Range button
    self._rangeButton = tk.Button(self.pointFrame, text="Set fit range", command=self._rangeButtonClicked)
    self._rangeButton.pack(side=tk.TOP, fill=tk.X)
    
    # Remove button
    self.removeButton = tk.Button(master=self.pointFrame, text="Remove component", \
                                  command=self._removeButtonClicked)
    self.removeButton.pack(side=tk.TOP, fill=tk.X)

    # a tk.DrawingArea
    self.canvas.get_tk_widget().pack()
    self.cid = self.f.canvas.mpl_connect('button_press_event', self._mouseButtonClicked)
    self.mwe = self.f.canvas.mpl_connect('scroll_event', self._mouseWheel)
    
    self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotFrame)
    self.toolbar.update()
    self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _quit():
      # stops main loop
      self.root.quit()
      # this is necessary on Windows to prevent
      # Fatal Python Error: PyEval_RestoreThread: NULL tstate
      self.root.destroy()

    self.quitButton = tk.Button(master=self.pointFrame, text='Quit', command=_quit)
    self.quitButton.pack(side=tk.BOTTOM, fill=tk.X)

  def _commonSigClicked(self):
    """
      Relate/untie the width parameter for Gaussian fit.
      
      By default, all widths are coupled to "sig1".
    """
    if self._compoCount < 2:
      return
    if self._commonSig.get() == 1:
      # Relate the widths
      for i in smo.range(1, self._compoCount+1):
        if i == 1:
          continue
        self._model.thaw("sig"+str(i))
        self._model.relate("sig"+str(i), ["sig"+str(1)])
    else:
      # Untie the widths
      after = "sig1" in self._model.freeParamNames()
      for i in smo.range(1, self._compoCount+1):
        if i == 1:
          continue
        self._model.untie("sig"+str(i), forceFree=after)

  def _mouseWheel(self, event):
    """
      React to mouse wheel
    """
    if self._compoCount < 1:
      return
    if self._parSelect.get() == 0:
      name = "off"
    else:
      name = self._compPars[self._parSelect.get()-1]
      name += str(self._activeComponent)
    
    if self._multAddSelect.get() == 1:
      # Multiplication
      try:
        if self._multNumberEntry.cget("bg") == "red":
          self._multNumberEntry.config(background=self.root.cget("bg"))
        delta = float(self._multNumber.get())
      except:
        self._multNumberEntry.config(background="red")
        return
      self._model[name] = self._model[name] * (1.0 + event.step*delta/100.0)
    else:
      # Add
      try:
        if self._addNumberEntry.cget("bg") == "red":
          self._addNumberEntry.config(background=self.root.cget("bg"))
        delta = float(self._addNumber.get())
      except:
        self._addNumberEntry.config(background="red")
        delta = self._config["defaultWheelAdd"]
      self._model[name] = self._model[name] + event.step*delta

    self._updateGUI()
    
  def _downAC(self):
    """
      Decrease no. of active component
    """
    if self._activeComponent > 1:
      self._activeComponent -= 1
    self._acstring.set(str(self._activeComponent))
    self._updateGUI()

  def _upAC(self):
    """
      Increase no. of active component
    """
    if self._activeComponent < self._compoCount:
      self._activeComponent += 1
    self._acstring.set(str(self._activeComponent))
    self._updateGUI()

  def _updateParValues(self):
    """
      Update GUI state of parameter values
    """
    if self._compoCount < 1:
      return
    self._parValues[0].set("%9.4f" % self._model["off"])
    for i, p in enumerate(self._compPars, 1):
      self._parValues[i].set("%9.4f" % (self._model[p+str(self._activeComponent)]))

  def _updateFreeCheckBox(self):
    """
      Update GUI state of free/frozen checkboxes
    """
    if self._compoCount < 1:
      return
    self._parFreeVars[0].set("off" in self._model.freeParamNames())
    for i, p in enumerate(self._compPars, 1):
      if (self._mode == "gauss") and (self._commonSig.get() == 1) and (p.find("sig") != -1):
        # Gaussian width is coupled so "sig1" is the relevant parameter.
        # Show all sigs as free of they are coupled and sig1 is free.
        self._parFreeVars[i].set("sig1" in  self._model.freeParamNames())
        continue
      self._parFreeVars[i].set((p+str(self._activeComponent)) in self._model.freeParamNames())

  def _freeCheckboxClicked(self):
    """
      Transfer GUI information to model
    """
    pars = ["off"]
    pars.extend([x+str(self._activeComponent) for x in self._compPars])
    
    if (self._mode == "gauss") and (self._commonSig.get() == 1):
      # If Gaussian widths are coupled, sig1 is the only relevant
      # parameter.
      for i, p in enumerate(pars):
        if p.find("sig") != -1:
          pars[i] = "sig1"
    
    free = []
    frozen = []
    for i, p in enumerate(pars):
      if self._parFreeVars[i].get() > 0:
        free.append(p)
      else:
        frozen.append(p)

    self._model.thaw(free)
    self._model.freeze(frozen)
  
  def _updateGUI(self):
    """
      Synchronize GUI and internal state.
    """
    self._updateModelPlot()
    self._updateView()
    self._updateParValues()
    self._updateFreeCheckBox()
    self._updateFreeCheckBox()
      
  def _fitModel(self):
    """
      Fit the model.
    """
    if self._yerr is None:
      self._model.fit(self._x[self._rangeIndices], self._y[self._rangeIndices])
    else:
      self._model.fit(self._x[self._rangeIndices], self._y[self._rangeIndices], \
                      yerr=self._yerr[self._rangeIndices])
    self._updateGUI()
    
  def _updateModelPlot(self):
    """
      Update the plot describing the current model.
    """
    # Remove existing lines
    if hasattr(self, "_modelLines"):
      for l in self._modelLines:
        self.a.lines.remove(l)
    
    self._modelLines = []
    model = self._model.evaluate(self._x)
    self.a.plot(self._x, model, self._config["modelLS"])
    self._modelLines.append(self.a.lines[-1])
    
    if self._config["plotIndComps"]:
      # Save parameters
      spars = self._model.parameters()
      # Set all areas to 0
      for i in smo.range(self._compoCount):
        n = str(i+1)
        self._model["A"+n] = 0.0
      # Add individual lines
      for i in smo.range(self._compoCount):
        n = str(i+1)
        self._model["A"+n] = spars["A"+n]
        imod = self._model.evaluate(self._x)
        # Assign line style
        if (i+1) == self._activeComponent:
          lst = self._config["activeCompLS"]
        else:
          lst = self._config["inactiveCompLS"]
        self.a.plot(self._x, imod, lst)
        self._modelLines.append(self.a.lines[-1])
        # Switch off the component again
        self._model["A"+n] = 0.0
      # Restore the old state of parameters
      self._model.assignValues(spars)
    
  def _cleanPointCache(self):
    """
      Empty point cache and remove lines from plot.
    """
    for p in self._pointCache:
      # Remove from plot
      self.a.lines.remove(p[1])
    self._pointCache = []
  
  def _clearRangeCache(self):
    """
      Empty range cache and remove any lines from plot.
    """
    for e in self._rangeCache:
      self.a.lines.remove(e[1])
    self._rangeCache = []
  
  def _rangeButtonClicked(self):
    """
      Initiate input of fit range.
    """
    if not hasattr(self, "_rangeCache"):
      self._rangeCache = []
    
    if (self._inputMode == "range") or (len(self._rangeCache) == 2):
      self._clearRangeCache()
      self._updateGUI()
    self._inputMode = "range"
  
  def _newCompAddPoint(self, event):
    """
      Add a point to point cache or create component if point set is complete.
    """
    p = Point(event)
    # Plot point on canvas
    self.a.plot([p.xdata], [p.ydata], 'r+', markersize=10)
    self._pointCache.append([p, self.a.lines[-1]])
    
    if len(self._pointCache) < 3:
      # More points have to be added 
      self._updateView()
      return
    
    # Three points have been added
    fwhm = np.abs(self._pointCache[0][0].xdata - self._pointCache[2][0].xdata)
    depth = self._model["off"] - self._pointCache[1][0].ydata
    # Estimate std (Gaussian)
    sigma = fwhm/2.35482
    # Estimate area of a Gaussian with given depth
    agauss = sigma*np.sqrt(2.*np.pi)*depth
    
    # Add further component
    self._compoCount += 1
    if self._mode == "voigt":
      nmodel = fuf.MultiVoigt1d(self._compoCount)
    elif self._mode == "gauss":
      nmodel = fuf.MultiGauss1d(self._compoCount)
    # Assign parameters from old model
    nmodel.assignValues(self._model.parameters())
    # Transfer free parameters 
    nmodel.thaw(self._model.freeParamNames())
    self._model = nmodel
    
    if self._mode == "voigt":
      c0 = 2.0056
      c1 = 1.0593
      # fl/fg
      phi = self._config["default_alad"]
      fg = fwhm/(1.0 - c0*c1 + np.sqrt(phi**2 + 2.*c1*phi + c0**2*c1**2))
      fl = phi*fg
      n = str(self._compoCount)
      self._model["A"+n] = -agauss
      self._model["al"+n] = fl/2.0
      self._model["ad"+n] = fg/(2.*(2.*np.log(2.0)))
      self._model["mu"+n] = self._pointCache[1][0].xdata
      # Thaw default parameters
      self._model.thaw([x+n for x in self._config["defaultFreeVoigt"]])
    elif self._mode == "gauss":
      n = str(self._compoCount)
      self._model["A"+n] = -agauss
      self._model["sig"+n] = sigma
      self._model["mu"+n] = self._pointCache[1][0].xdata
      # Thaw default parameters
      self._model.thaw([x+n for x in self._config["defaultFreeGauss"]])
      if (self._commonSig.get() == 1) and (self._compoCount > 1):
        # Common widths for Gaussians requested
        self._model.thaw("sig"+n)
        self._model.relate("sig"+n, ["sig"+str(1)])
    
    self._cleanPointCache()
    self._updateGUI()
    
  def _addRangePoint(self, event):
    """
      Define fit range and finalize of complete.
    """
    p = Point(event)
    # Plot point on canvas
    self.a.plot([p.xdata]*2, [self._y.min(), self._y.max()], self._config["rangeLS"])
    self._rangeCache.append([p, self.a.lines[-1]])
    self._updateGUI()
    
    if len(self._rangeCache) < 2:
      # Another point has to be added
      return
    
    # There are two points. Therefore, the range can be set.
    self._range = [self._rangeCache[0][0].xdata, self._rangeCache[1][0].xdata]
    self._rangeIndices = np.where((self._x >= self._range[0]) & (self._x <= self._range[1]))[0]
    # Reset mode
    self._inputMode = "newComp"
    
  def _mouseButtonClicked(self, event):
    """
      Called on click of mouse button.
    """
    # Accept only middle button
    if event.button != 2: return
    # Check whether click occurred on plot axes.
    if event.inaxes is not self.a: return
    
    if self._inputMode == "newComp":
      self._newCompAddPoint(event)
    elif self._inputMode == "range":
      self._addRangePoint(event)
    
  def _removeButtonClicked(self):
    """
      Remove the active component.
    """
    if self._compoCount < 2:
      return
    
    # Create new model
    if self._mode == "voigt":
      nm = fuf.MultiVoigt1d(self._compoCount-1)
    else:
      nm = fuf.MultiGauss1d(self._compoCount-1)
    
    # The offset has to be treated separately
    nm["off"] = self._model["off"]
    if "off" in self._model.freeParamNames():
      nm.thaw("off")
    
    # t counts the "target" component
    t = 0
    for i in smo.range(1, self._compoCount+1):
      if i == self._activeComponent:
        # Disregard the currently active component
        continue
      n = str(i)
      t += 1
      for p in self._compPars:
        nm[p + str(t)] = self._model[p + n]
        if (p+n) in self._model.freeParamNames():
          nm.thaw(p + str(t))

    self._model = nm
    self._compoCount -= 1
    self._downAC()
    self._updateGUI()
      
  def _updateView(self):
    """
      Redraw MPL canvas
    """
    self.f.canvas.draw()

  def interactiveFit(self):
    """
      Carry out the interactive fit
      
      Returns
      -------
      Parameters : dict
          Dictionary mapping parameter name to the value.
      Active parameters : dict
          Dictionary mapping active parameter names (not including
          any numbers) to value
      Active component : int
          Number of the 'active component'. Counting starts at one.
    """
    self.root.wm_title(self.windowTitle)
    self.canvas.show()
    tk.mainloop()
    # Prepare return value
    aps = {"off":self._model["off"]}
    n = str(self._activeComponent)
    for p in self._compPars:
      aps[p] = self._model[p+n]
    return self._model.parameters(), aps, self._activeComponent
    
    
