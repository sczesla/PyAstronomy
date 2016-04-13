import numpy as np
import scipy.interpolate as sci
import gzip

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyasl
from bisect import bisect

import six.moves.tkinter as tk
import six

from .pyaPicker import Point


class ContinuumInteractive:
  """
    GUI for interactive point selection.
    
    The `ContinuumInteractive` class provides a tool to normalize data
    interactively. The data are plotted and a number of points are defined
    using the mouse, which are connected by a spline to define a
    continuum estimate.
    
    Points are *selected* using the *middle button* of the mouse on the
    plot window. 
    
    Parameters
    ----------
    x, y : array
        The data to be normalized.
    config : dictionary, optional
        Information used for configuration:
          - specPlotStyle: Style used to plot the data (default 'b.--')
          - astyle: Style used to plot 'active' points (default 'ro')
          - istyle: Style used to plot 'inactive' points (default 'yp')
          - splineLineStyle: Style used to plot the continuum estimate (default 'r--')
          - normPlotRefLine: If True (default), a reference line is printed in the normalization window 
          - normLineStyle: Style used to plot the normalized data (default 'b.--')
          - normRefLineStyle: Style used to plot the reference line for the normalization (default 'k--')
          - sortPointListX: Determines whether point list is sorted in ascending order (default: True)
          - usegzip: If True (default), gzip is used to open and save files.
          - windowTitle: Title of the window
    
    Attributes
    ----------
    f : mpl Figure
        A Figure instance from matplotlib.
    a : mpl axis
        An Axis instance from matplotlib
    pointList : list of points
        A list containing the selected points in the form
        of `Point` instances.
    astyle : string
        Mpl style for plotting used for "active" point.
        Default is "ro".
    istyle : string
        Mpl style for plotting used for "inactive" point.
        Default is "yp".
  """
  
  def __init__(self, x, y, config=None):
    
    dconfig = {"specPlotStyle":"b.--", 
              "astyle":"ro",
              "istyle":"yp",
              "splineLineStyle":"r--",
              "normLineStyle":"b.--",
              "normPlotRefLine":True,
              "normRefLineStyle":"k--",
              "sortPointListX":True,
              "usegzip":True,
              "windowTitle":"PyA Continuum Interactive"}
    if config is not None:
      for k in six.iterkeys(dconfig):
        if k in config:
          dconfig[k] = config[k]

    self.config = dconfig
    
    self.windowTitle = dconfig["windowTitle"]
    self.f = Figure()
    self.a = self.f.add_subplot(111)
    
    # The normalized plot figure
    self.normf = Figure()
    self.norma = self.normf.add_subplot(111)
    self._normaLineRef = None
    
    # Save the data (spectrum)
    self._x = x.copy()
    self._y = y.copy()
    self.a.plot(x, y, dconfig["specPlotStyle"])
    
    # Active and Inactive plot style for points
    self.astyle = dconfig["astyle"]
    self.istyle = dconfig["istyle"]
    
    self.pointList = []

    self.root = tk.Tk()
    
    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.pack(fill=tk.BOTH, side=tk.LEFT, expand=True)
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.pointFrame = tk.Frame(self.root)
    self.pointFrame.pack(side=tk.LEFT, fill=tk.BOTH)
    
    self.listLabel = tk.Label(self.pointFrame, text="Selected points")
    self.listLabel.pack(side=tk.TOP)
    
    self.lb = tk.Listbox(self.pointFrame, height=15, selectmode=tk.SINGLE)
    self.lb.pack(side=tk.TOP)
    self.lb.bind("<<ListboxSelect>>", self._lbSelect)
    
    self.removeButton = tk.Button(master=self.pointFrame, text="Remove", \
                                  command=self._removeButtonClicked)
    self.removeButton.pack(side=tk.TOP, fill=tk.X)
    
    self.clearButton = tk.Button(master=self.pointFrame, text="Clear all", \
                                 command=self._clearAll)
    self.clearButton.pack(side=tk.TOP, fill=tk.X)

    # Says whether the normalized spectrum is currently shown 
    self._normalizedDataShown = False
    # The dummy label only holds some space
    self.dummyLabel1 = tk.Label(self.pointFrame, text="")
    self.dummyLabel1.pack(side=tk.TOP)
    self.splineBoxLabel = tk.Label(self.pointFrame, text="Spline kind")
    self.splineBoxLabel.pack(side=tk.TOP)
    # Define spline selection
    self.splineOptions = ["Linear", "Quadratic", "INTEP"]
    self.splineFuncs = [self._evalLinearSpline, self._evalQudraticSpline, self._evalIntepSpline]
    self.splineSelectBox = tk.Listbox(self.pointFrame, height=len(self.splineOptions), selectmode=tk.SINGLE, \
                                      exportselection=False)
    for item in self.splineOptions:
      self.splineSelectBox.insert(tk.END, item)
    self.splineSelectBox.pack(side=tk.TOP)
    self.splineSelectBox.bind("<<ListboxSelect>>", self._splineSelected)
    self.splineSelectBox.select_set(0)
    self._splineSelected(None)

    self.normButton = tk.Button(master=self.pointFrame, text="Normalize", \
                                  command=self._showNorm)
    self.normButton.pack(side=tk.TOP, fill=tk.X)

    # a tk.DrawingArea
    self.canvas.get_tk_widget().pack()
    self.cid = self.f.canvas.mpl_connect('button_press_event', self._mouseButtonClicked)
    
    self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotFrame)
    self.toolbar.update()
    self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _quit():
      if self._normalizedDataShown:
        self._quitWin()
      # stops main loop
      self.root.quit()
      # this is necessary on Windows to prevent
      # Fatal Python Error: PyEval_RestoreThread: NULL state
      self.root.destroy()

    self.quitButton = tk.Button(master=self.pointFrame, text='Quit', command=_quit)
    self.quitButton.pack(side=tk.BOTTOM, fill=tk.X)

    self.saveLoadFrame = tk.Frame(self.pointFrame)
    self.saveLoadFrame.pack(side=tk.BOTTOM, fill=tk.X)
    self.saveButton = tk.Button(master=self.saveLoadFrame, text="Save", command=self._saveToFile)
    self.saveButton.pack(side=tk.LEFT, fill=tk.X, expand=True)
    self.loadButton = tk.Button(master=self.saveLoadFrame, text="Load", command=self._loadFromFile)
    self.loadButton.pack(side=tk.RIGHT, fill=tk.X, expand=True)
    
    self.root.protocol("WM_DELETE_WINDOW", _quit)


  def _fileOpenMethod(self):
    """
      Get method for opening file.
      
      Determines whether gzip is used for file opening or not.
      
      Returns
      -------
      open : callable
          Method to open file.
    """
    if self.config["usegzip"]:
      return gzip.open
    return open

  def _pickleFileExtension(self):
    """
      Get default extension for pickle files.
      
      Returns
      -------
      Extension : string
          The default file extension.
    """
    if self.config["usegzip"]:
      return ".pickle.gz"
    else:
      return ".pickle"

  def _saveToFile(self, fn=None):
    """
      Save state to a pickle file.
    """
    import pickle
    pfe = self._pickleFileExtension()
    if fn is None:
      # Request a filename
      import six.moves.tkinter_tkfiledialog as tkFileDialog
      fn = tkFileDialog.asksaveasfilename(defaultextension=pfe, title="Save as pickle file", \
                                          filetypes=[("pickle files", "*"+pfe)])
    pickle.dump(self._getState(), self._fileOpenMethod()(fn, 'w'))

  def saveStateToPickleFile(self, fn):
    """
      Save the state (point selection etc.) to pickle file.
      
      Parameters
      ----------
      fn : string
          The filename.
    """
    self._saveToFile(fn=fn)

  def _loadFromFile(self):
    """
      Load state from a pickle file.
    """
    import pickle
    import six.moves.tkinter_tkfiledialog as tkFileDialog
    pfe = self._pickleFileExtension()
    fn = tkFileDialog.askopenfilename(defaultextension=pfe, title="Load from pickle file", \
                                      filetypes=[("pickle files", "*"+pfe)])
    state = pickle.load(self._fileOpenMethod()(fn))
    
    class Event:
      def __init__(self, xd, yd, inax):
        self.button = 2
        self.inaxes = inax
        self.xdata = xd
        self.ydata = yd
    
    for p in state["points"]:
      e = Event(p[0], p[1], self.a)
      self._mouseButtonClicked(e)
    
    self.splineSelectBox.select_clear(0, len(self.splineOptions))
    selected = self.splineOptions.index(state["splineKind"])
    self.splineSelectBox.select_set(selected)
    self._splineSelected(None)

  
  def _updateNormalizationPlot(self):
    """
      Replot the normalized spectrum.
    """
    if not self._normalizedDataShown:
      return
    
    normFlux = self._y / self._currentSpline
    
    if not self._normaLineRef is None:
      self.norma.lines.pop(self.norma.lines.index(self._normaLineRef))
    else:
      # First plot, add the reference line if requested
      if self.config["normPlotRefLine"]:
        self.norma.plot([self._x[0], self._x[-1]], [1.0, 1.0], self.config["normRefLineStyle"])
      
    self.norma.plot(self._x, normFlux, self.config["normLineStyle"])
    self._normaLineRef = self.norma.lines[-1]
    self.normCanvas.draw()
    
  def _quitWin(self):
    self._normalizedDataShown = False
    if not self._normaLineRef is None:
      self.norma.lines.pop(self.norma.lines.index(self._normaLineRef))
    self._normaLineRef = None
    self._normwin.destroy()
  
  def _showNorm(self):
    """
      Shows normalized data in a separate window.
    """
    
    if self._normalizedDataShown:
      return
    
    self._normwin = tk.Tk()
    self._normwin.wm_title("Normalized spectrum")
    self._normwin.protocol("WM_DELETE_WINDOW", self._quitWin)
      
    # A frame containing the mpl plot
    self.normFrame = tk.Frame(master=self._normwin)
    self.normFrame.pack(fill=tk.BOTH, side=tk.LEFT, expand=True)
    self.normCanvas = FigureCanvasTkAgg(self.normf, master=self.normFrame)

    # a tk.DrawingArea
    self.normCanvas.get_tk_widget().pack()
    self.normToolbar = NavigationToolbar2TkAgg(self.normCanvas, self.normFrame)
    self.normToolbar.update()
    self.normCanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
      
    closeButton = tk.Button(master=self._normwin, text='Close', command=self._quitWin)
    closeButton.pack(side=tk.BOTTOM)
    
    self._normalizedDataShown = True
    self._updateNormalizationPlot()
  
  def _pointListToArray(self):
    """
      Convert pointList data values into x,y arrays
      
      Returns
      -------
      x, y : array
          The x/y-data for the selected points. 
    """
    x, y = np.zeros(len(self.pointList)), np.zeros(len(self.pointList))
    for i, p in enumerate(self.pointList):
      x[i], y[i] = p.xdata, p.ydata
    indi = np.argsort(x)
    return x[indi], y[indi] 
  
  def _evalLinearSpline(self):
    """
      Evaluates linear spline.
      
      Returns
      -------
      Spline : array
          The continuum estimate.
    """
    x, y = self._pointListToArray()
    if len(x) < 2:
      return self._x * np.nan
    f = sci.interp1d(x, y, bounds_error=False, fill_value=np.nan)
    return f(self._x)
  
  def _evalQudraticSpline(self):
    """
      Evaluates quadratic spline.
      
      Returns
      -------
      Spline : array
          The continuum estimate.
    """
    x, y = self._pointListToArray()
    if len(x) < 3:
      return self._x * np.nan
    f = sci.interp1d(x, y, bounds_error=False, fill_value=np.nan, kind="quadratic")
    return f(self._x)
    
  def _evalIntepSpline(self):
    """
      Evaluates INTEP spline.
      
      Returns
      -------
      Spline : array
          The continuum estimate.
    """
    x, y = self._pointListToArray()
    if len(x) < 2:
      return self._x * np.nan
    return pyasl.intep(x, y, self._x, boundsError=False, fillValue=np.nan)
    
  def _splineSelected(self, event):
    """
      Called when the spline selection box is clicked.
    """
    selected = int(self.splineSelectBox.curselection()[0])
    self._evalSpline = self.splineFuncs[selected]
    self._evalSpline()
    self._updateSplineLine()
    self._updateView()
    # Save the name of the current spline selection
    self._currentSplineKind = self.splineOptions[selected]
    
  def _plotPointAs(self, state, point=None, lbString=None):
    """
      Plot a point in active/inactive color
      
      Parameters
      ----------
      state : string, {"active", "inactive"}
          Which color to use
      point : Point, optional
          The point information as an instance of the Point class.
      lbString : string
          The point identifier used in the listbox.
    """
    if (point is None) == (lbString is None):
      raise(PE.PyAValError("Either `point` or `lbString` have to be specified."))
    if point is not None:
      pli, lli = self._searchPoint(point.lbIdent)
    else:
      pli, lli = self._searchPoint(lbString)
    # Remove 'old' point from plot
    self.a.lines.pop(lli)
    if state == "active":
      style = self.astyle
    elif state == "inactive":
      style = self.istyle
    # Add new point in specified color
    self.a.plot([self.pointList[pli].xdata], [self.pointList[pli].ydata], style)
    self.pointList[pli].mplLine = self.a.lines[-1]
      
  def _searchPoint(self, lbString):
    """
      Search point specified by listbox string. 
      
      Parameters
      ----------
      lbString : string
          The string used as identifier in the listbox.
      
      Returns
      -------
      Point list index : int
          The index of the point in the `pointlist`.
      Line index : int
          The index of the corresponding line in the
          lines list attribute of the axis instance. 
    """
    for ip, p in enumerate(self.pointList):
      if p.lbIdent == lbString: break
    for i, l in enumerate(self.a.lines):
      if l is p.mplLine:
        break
    return ip, i
    
  def _lbSelect(self, event):
    """
      React on change of point-listbox selection.
      
      This function determines the selected item
      and plots the point in "active" style.
    """
    # Get 'old' active item (changed afterward)
    oa = self.lb.get(tk.ACTIVE)
    # Get active item (string)
    ai = self.lb.get(int(self.lb.curselection()[0]))
    self._plotPointAs("inactive", lbString=oa)
    self._plotPointAs("active", lbString=ai)
    # Update view
    self._updateView()
    
  def _updateSplineLine(self):
    """
      Draws a new line representing the spline.
    """
    if not hasattr(self, "_splineLineRef"):
      self._splineLineRef = None
    if not self._splineLineRef is None:
      # Remove old line
      self.a.lines.pop(self.a.lines.index(self._splineLineRef))
    self._currentSpline = self._evalSpline()
    self.a.plot(self._x, self._currentSpline, self.config["splineLineStyle"])
    self._splineLineRef = self.a.lines[-1]
    if self._normalizedDataShown:
      self._updateNormalizationPlot()
    
  def _mouseButtonClicked(self, event):
    """
      Called on click of mouse button.
      
      If the middle button has been clicked, a
      point is added to the selection.
    """
    # Accept only middle button
    if event.button != 2: return
    # Check whether click occurred on plot axes.
    if event.inaxes is not self.a: return
    
    p = Point(event)
    self.a.plot([p.xdata], [p.ydata], self.istyle)
    p.mplLine = self.a.lines[-1]
    
    index = tk.END
    if self.config["sortPointListX"]:
      # Insert point into the list at appropriate position
      # to keep x-values sorted
      x, y = self._pointListToArray()
      index = bisect(x, p.xdata)
    
    # Insert the point's content into the list box
    self.lb.insert(index, p.asStr())
    # Save the "list box identifier"
    p.lbIdent = p.asStr()
    # Save the point
    self.pointList.append(p)
    
    self._updateSplineLine()
    self._updateView()

  def _removeButtonClicked(self):
    """
      Remove a point.
    """
    # Check whether anything was selected
    sel = self.lb.curselection()
    if len(sel) == 0: return
    # Find indices for the point 
    el = self.lb.get(int(sel[0]))
    pli, lli = self._searchPoint(el)
    # Remove it
    self.a.lines.pop(lli)
    self.pointList.pop(pli)
    self.lb.delete(sel[0])
    self._updateSplineLine()
    self._updateView()

  def _clearAll(self):
    """
      Remove all previously selected points.
    """
    np = self.lb.size()
    for _ in six.moves.range(np):
      self.lb.select_clear(0, self.lb.size())
      self.lb.select_set(0)
      self._removeButtonClicked()

  def _updateView(self):
    """
      Redraw MPL canvas
    """
    self.f.canvas.draw()

  def plot(self, *args, **kwargs):
    """
      Plot on interactive canvas.
      
      Accepts all arguments and keywords also accepted by
      matplotlib's `plot` method.
    """
    self.a.plot(*args, **kwargs)

  def plotNorm(self, *args, **kwargs):
    """
      Plot on normalization window.
      
      Accepts all arguments and keywords also accepted by
      matplotlib's `plot` method.
    """
    self.norma.plot(*args, **kwargs)

  def _getState(self):
    """
      Collect current state in a dictionary.
      
      Returns
      -------
      State : dictionary
          The following keys are defined:
            - points: A list of two-float tuples holding the
                      x,y location of the selected points.
            - continuum : Array holding the continuum estimate
                          at the given x-values.
            - splineKind : A string specifying the selected
                           spline option.
            - normalizedData : An array holding the normalized
                               data.
    """
    result = {}
    result["points"] = []
    for p in self.pointList:
      result["points"].append((p.xdata, p.ydata))
    result["continuum"] = self._currentSpline.copy()
    result["normalizedData"] = self._y / result["continuum"]
    result["splineKind"] = self._currentSplineKind
    return result

  def findContinuum(self):
    """
      Interactively find the continuum estimate.

      Returns
      -------
      State : dictionary
          The following keys are defined:
            - points: A list of two-float tuples holding the
                      x,y location of the selected points.
            - continuum : Array holding the continuum estimate
                          at the given x-values.
            - splineKind : A string specifying the selected
                           spline option.
            - normalizedData : An array holding the normalized
                               data.
    """
    self.root.wm_title(self.windowTitle)
    self.canvas.show()
    tk.mainloop()
    return self._getState()
    

