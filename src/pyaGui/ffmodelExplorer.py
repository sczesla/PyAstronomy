# -*- coding: utf-8 -*-
from __future__ import print_function, division
import matplotlib

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from PyAstronomy.pyaC import pyaErrors as PE
import numpy
import six.moves.tkinter as tk
import six.moves.tkinter_messagebox as tkMessageBox
import platform

class FFModelPlotFit:
  """
    Plotting and fitting for the model.
    
    Parameters
    ----------
    x : array
        The abscissa values.
    y : array
        The ordinate values.
    yerr : array, optional
        The error on the points.
    withResiduals : boolean, optional
        If True, the plot will also show the
        residuals. If `yerr` is given, residuals
        will be plotted in units of standard
        deviations. Default is False.
  """
  
  def __init__(self, x, y, yerr=None, withResiduals=False):
    self.x = x
    self.y = y
    self.yerr = yerr
    self.a = None
    self.withResiduals = withResiduals
    self.modelLine = None
    self.residualLine = None
    self.dataPlot = False
    self.fitLims = None # range for fit
    self.fitIdx = None # indices for fit
    self.chi2 = None
    
  def _setUpAxes(self, f):
    """
      Create axes for plotting
    """
    if not self.withResiduals:
      self.a = f.add_subplot(111)
    else:
      gs = matplotlib.gridspec.GridSpec(3, 1)
      self.a = f.add_subplot(gs[0:2,::])
      self.ar = f.add_subplot(gs[2,::], sharex=self.a)
  
  def plot(self, f, odf):
    """
      Plot the model.
      
      Parameters
      ----------
      f : matplotlib.Figure
          The figure on which to plot.
      odf : funcFit model class
          The model which is explored.
    """
    # Record whether this is the first call to adjust
    # the autoscale parameter
    firstCall = False
    if self.a is None:
      self._setUpAxes(f)
      firstCall = True
    
    # Calculate the model
    model = odf.evaluate(self.x)
    
    if not self.dataPlot:
      # Plot the data
      self.dataPlot = True
      if self.yerr is None:
        self.a.plot(self.x, self.y, 'bp')
        self.a.plot(self.x, self.y, ms='',ls='-',color='b',alpha=0.2)
      else:
        self.a.errorbar(self.x, self.y, yerr=self.yerr, fmt='b.')
    
    if self.modelLine is not None:
      # Remove the last model to show the latest one
      self.modelLine[0].remove()
    # Plot the model
    self.modelLine = self.a.plot(self.x, model, 'r--')
    
    # Derive chi2 value    
    self.chi2 = self.get_chi2(odf)
    
    if self.withResiduals:
      if self.residualLine is not None:
        self.ar.lines.pop(-1)
      res = self.y - model
      if self.yerr is not None:
        res /= self.yerr
      self.residualLine = self.ar.plot(self.x, res, 'bp')
    if firstCall:
      self.a.autoscale(False)
  
  def get_chi2(self, odf, lims = None):
    """
      Calculate current chi2 value
    """
    if lims:
        self.fitLims = lims

    r = 0.0
    model = odf.evaluate(self.x)
    if self.fitLims:
        #print "Using lims: ",self.fitLims
        self.fitIdx = numpy.where(numpy.logical_and(self.x >= self.fitLims[0], self.x <= self.fitLims[1]))[0]
    else:
        self.fitIdx = numpy.arange(len(self.x))
    if self.yerr is not None:
        r = numpy.sum((self.y[self.fitIdx]-model[self.fitIdx])**2 / self.yerr[self.fitIdx]**2)
    else:
        r = numpy.sum((self.y[self.fitIdx]-model[self.fitIdx])**2)
    return r

  def fit(self, odf, lims = None, verbose=1):
    """
      Fit the model.
      
      Parameters
      ----------
      odf : funcFit model class
          The model to be fitted.
    """
    
    if lims:
        if False: # if old fit ranges -> use old indices... lims == self.fitLims and self.fitIdx:
            odf.fit(self.x[self.fitIdx], self.y[self.fitIdx], yerr=self.yerr[self.fitIdx])
        else:
            self.fitLims = lims
            self.fitIdx = numpy.where(numpy.logical_and(self.x >= self.fitLims[0], self.x <= self.fitLims[1]))[0]
            if verbose > 1:
              print("fit index: ",self.fitIdx)
              print("         x ",self.x[self.fitIdx])
              print("         y ", self.y[self.fitIdx])
              if self.yerr: print("    yerr ",self.yerr[self.fitIdx])
              print("     model ",odf.evaluate(self.x[self.fitIdx]))
            if self.yerr is not None:
              odf.fit(self.x[self.fitIdx], self.y[self.fitIdx], yerr=self.yerr[self.fitIdx])
            else:
              odf.fit(self.x[self.fitIdx], self.y[self.fitIdx])
    else:
        if self.yerr is not None:
          odf.fit(self.x, self.y, yerr=self.yerr)
        else:
          odf.fit(self.x, self.y)


    
class SetToDialog(tk.Toplevel):
  
  def __init__(self, parent, oldVal, pname):
    tk.Toplevel.__init__(self, parent)
    self.transient(parent)
    self.parent = parent
    self.newVal = None
    
    self.inputFrame = tk.Frame(self)
    self.parLabel = tk.Label(self.inputFrame, text="Value for parameter: " + str(pname))
    self.parLabel.pack()
    self.inputVal = tk.StringVar()
    self.inputVal.set("% g" % oldVal)
    self.input = tk.Entry(self.inputFrame, textvariable=self.inputVal)
    self.input.pack()
    self.inputFrame.pack(fill=tk.X)
    
    self.buttonFrame = tk.Frame(self)
    self.okButton = tk.Button(self, text="OK", command=self._okClicked)
    self.caButton = tk.Button(self, text="Cancel", command=self._cancelClicked)
    self.okButton.pack(side=tk.RIGHT)
    self.caButton.pack(side=tk.LEFT)
    
    # Treat return as OK
    self.bind("<Return>", self._okClicked)
    # Treat close as cancel
    self.protocol("WM_DELETE_WINDOW", self._cancelClicked)
    # For a modal dialog
    self.grab_set()
    
    # Keyboard focus to input entry
    self.input.focus_set()
    
    self.wait_window(self)
  
  def _okClicked(self, *args):
    try:
      mf = float(self.inputVal.get())
    except ValueError:
      tkMessageBox.showwarning("Invalid float", "Cannot convert " + self.modEntry.get() + " to float." + \
                               " Make it a valid value to proceed.")
      return
    self.newVal = mf
    self._cancelClicked()
  
  def _cancelClicked(self, *args):
    self.parent.focus_set()
    self.destroy()


class ShowParameterSummary(tk.Toplevel):
  
  def __init__(self, parent, odf):
    tk.Toplevel.__init__(self, parent)
    self.wm_title("Parameter Summary")
    self.parent = parent
    self.newVal = None
    
    self.showFrame = tk.Frame(self)
    self.showText = tk.Text(self.showFrame)
    self.showText.config(background='white')
    self.showText.pack(fill=tk.BOTH, expand=True)
    self.showFrame.pack(fill=tk.BOTH, expand=True)
    
    self.protocol("WM_DELETE_WINDOW", self._windowClosed)
    
  def updateInfo(self, odf):
    """
      Update the parameter information.
    """
    self.showText.config(state=tk.NORMAL)
    try:
      self.showText.delete("1.0", tk.END)
    except tk.TclError:
      # Ignore the absence of text 
      pass
    lines = odf.parameterSummary(toScreen=False)
    self.showText.config(width=max([len(x) for x in lines]), \
                          height=len(lines)+1)

    t = "".join([x + "\n" for x in lines])
    self.showText.insert(tk.END, t)
    self.showText.config(state=tk.DISABLED)

  def _windowClosed(self):
    """
      The window is closed.
    """
    # Tell the parent that the window does no longer exist
    self.parent.parSumWin = None
    self.destroy()


class ValueSetSummary(tk.Toplevel):
  """
    Inspired from ShowParameterSummary
  """
  def __init__(self, parent, odf):
    tk.Toplevel.__init__(self, parent)
    self.wm_title("Value set code")
    self.parent = parent
    self.newVal = None
    self.odf = odf
    self.showFrame = tk.Frame(self)
    
    self.textLabel = tk.Label(self.showFrame, text="Model name: ")
    self.textLabel.grid(row=0,column=0)
    self.modelName = tk.StringVar()
    self.modelName.set("model")
    self.modelEntry = tk.Entry(self.showFrame, textvariable=self.modelName, width=8)
    self.modelEntry.bind('<Return>', self._updateModelName)
    self.modelEntry.grid(row=0,column=1)
    
    self.showText = tk.Text(self.showFrame)
    self.showText.config(background='white')
    self.showText.grid(row=1,column=0,columnspan=2) #.pack(expand=True)
    
    self.clipboardButton = tk.Button(self.showFrame, text="Copy to Clipboard", command=self._copyToClipboard)
    self.clipboardButton.grid(row=2, column=0, columnspan=2)
    
    self.helpLabel = tk.Label(self.showFrame, text="You may use the above code for pasting...")
    self.helpLabel.grid(row=3, column=0, columnspan=2)
    
    self.showFrame.pack()
    
    self.protocol("WM_DELETE_WINDOW", self._windowClosed)

  def _updateModelName(self, *args):
      self.updateInfo(self.odf)

  def updateInfo(self, odf):
    """
      Update the parameter information.
    """
    self.showText.config(state=tk.NORMAL)
    try:
      self.showText.delete("1.0", tk.END)
    except tk.TclError:
      # Ignore the absence of text 
      pass
    
    lines = []
    for p in odf.availableParameters():
        lines.append(self.modelName.get()+"[\""+p+"\"] = "+str(odf[p]))
    
    self.showText.config(width=max([len(x) for x in lines]), \
                          height=len(lines)+1)
    t = "".join([x + "\n" for x in lines])
    self.showText.insert(tk.END, t)
    self.showText.config(state=tk.DISABLED)

  def _copyToClipboard(self):
    self.clipboard_clear()
    text = self.showText.get("1.0",tk.END)
    self.clipboard_append(text)
  def _windowClosed(self):
    """
      The window is closed.
    """
    # Tell the parent that the window does no longer exist
    self.parent.parSumWin = None
    self.destroy()


class FFModelExplorerList:
  """
    Adapt model interactively.

    Parameters
    ----------
    odf : Instance of OneDFit
        The model to be adapted.
    plotter : Instance of FFModelPlotFit or custom
        Class instance managing the actual plotting.
  """
  
  def __init__(self, odf, plotter):
    self.f = Figure()
    
    # Save reference to the model
    self.odf = odf
    # Save reference to the plotter
    self.plotter = plotter
    
    self.leftMask = None
    self.rightMask = None
    self.activeLine = None

    self.root = tk.Tk()
    self.root.wm_title("PyA Model Explorer")
    # Make the widgets expand/shrink as window size changes
    self.root.columnconfigure(0, weight=1)
    self.root.rowconfigure(0, weight=1)
    
    # Bind the mouse wheel
    if platform.system() == "Linux":
        self.root.bind("<Button-4>", self._mouseWheel)
        self.root.bind("<Button-5>", self._mouseWheel)    
    elif platform.system() == "Darwin":
        self.root.bind("<MouseWheel>", self._onWheel) # for OS X

    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.pack(fill=tk.BOTH, side=tk.LEFT, expand=True)
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.controlFrame = tk.Frame(self.root)
    self.controlFrame.pack(side=tk.RIGHT, expand=False, fill=tk.BOTH)
    self.selectedPar = tk.StringVar(self.controlFrame)

    # Get parameters of model
    ps = list(self.odf.parameters().keys())
    # Set default modification properties
    # Saves these properties for all parameters
    self.modProps = {}
    for p in ps:
      self.modProps[p] = {"modus":"mul", "modValMul":1.02, "modValAdd":0.01}
    
    # Frame containing all parameters
    self.parameterFrame = tk.Frame(self.controlFrame, height=2, bd=1, relief=tk.SUNKEN)    
    
    # Dictionaries holding the specific information for each parameter
    self.singleParameterFrames = {}
    self.singleParameterEntry = {}
    self.singleParameterVar = {}
    self.singleParameterFree = {} # knows whether the parameter is free (value=True) or frozen (=False)
        
    # Closures around the functions adapting thaw/freeze    
    def frozenChanged(k):
        def change():
            if self.singleParameterFree[k].get() == True:
                self.odf.thaw(k)
            else:
                self.odf.freeze(k)
            self._updateDof()     
        return change
    
    # define what happens when a value is set to the entered parameter (closures)
    def hitReturn(k):
      def change(*args):
        self.selectedPar.set(k)
        self.odf[k] = float(self.singleParameterVar[k].get())
        self._parameterValueChanged()
      return change  
    
    # defines what happens when a parameter's value is changed, but not set yet, i.e., not current (closures)
    def parameterValueChanged(k):
        def valueChanged(*args):
            pp,ll = str(self.singleParameterVar[k].get()).find("."), len(str(self.singleParameterVar[k].get()))
            if round(self.odf[k],ll-pp-1) != round(self.singleParameterVar[k].get(), ll-pp-1):
                self.singleParameterEntry[k].configure(bg = "red")
            else:
                self.singleParameterEntry[k].configure(bg = "white")
        return valueChanged    
    
    
    # Create an entry for each parameter
      ## Create a scrollable region
      ## Check maximum number of characters for parameter names:
    maxParamLen = 0
    for k in sorted(self.odf.parameters().keys()):
      if len(k) > maxParamLen:
        maxParamLen = len(k)
      ## Create an entry for each parameter	
    for k in sorted(self.odf.parameters().keys()):
      x = tk.Frame(self.parameterFrame, height=2, bd=2, relief=tk.SUNKEN,pady=2)
      self.singleParameterFrames[k] = x 
      y0 = tk.Radiobutton(x, text=k,variable=self.selectedPar,value=k, width=maxParamLen+1, indicatoron=0)
      y0.pack(side=tk.LEFT)
      #y1 = tk.StringVar()
      y1 = tk.DoubleVar()
      y1.set(self.odf[k])
      y2 = tk.Entry(x, textvariable=y1, width=8)
      y2.bind('<Return>', hitReturn(k))
      self.singleParameterVar[k] = y1
      self.singleParameterVar[k].trace('w', parameterValueChanged(k))
      self.singleParameterEntry[k]= y2
      y2.pack(side=tk.LEFT)
      self.singleParameterFrames[k].pack()
      modModus = tk.StringVar()
      modModus.set(self.modProps[k]["modus"])
      self.singleParameterFree[k] = tk.BooleanVar()
      self.singleParameterFree[k].set(k in self.odf.freeParamNames())
      y3 = tk.Radiobutton(x, text="Thawed", value=True,variable=self.singleParameterFree[k], command=frozenChanged(k))
      y4 = tk.Radiobutton(x, text="Frozen", value=False,variable=self.singleParameterFree[k], command=frozenChanged(k))
      y3.pack(side=tk.RIGHT)
      y4.pack(side=tk.RIGHT)
      self.parameterFrame.pack(fill=tk.X, expand=False)
      
    # Set of the menu to select the current parameter
    self.selectedPar.set(ps[0])
    
    
    # Chi2 frame:
    self.chi2Frame = tk.Frame(self.controlFrame, borderwidth=10)#, relief=tk.SUNKEN)
    self.chi2value = tk.DoubleVar()
    self.chi2value.set(self.plotter.chi2)
    self.dofValue = tk.IntVar()
    #self.dofValue.set(len(self.plotter.fitIdx))
    #self.dofValue.set(100)
    self.dofValue.set(len(self.plotter.x) - len(self.odf.freeParameters())-1)
    self.chi2Label = tk.Label(self.chi2Frame,text="Chi2: ")
    self.dofLabel = tk.Label(self.chi2Frame,text="dof: ")
    self.chi2Entry = tk.Entry(self.chi2Frame, textvariable=self.chi2value, bd=2, width=10)
    self.dofEntry = tk.Entry(self.chi2Frame, textvariable=self.dofValue, bd=2, width=10)
    self.chi2Label.pack(side=tk.LEFT)
    self.chi2Entry.pack(side=tk.LEFT)
    self.dofLabel.pack(side=tk.LEFT)
    self.dofEntry.pack(side=tk.LEFT)
    self.chi2Frame.pack()
  
    
    # Frame to bundle mouse-wheel inputs
    self.mouseWheelFrame = tk.Frame(self.controlFrame, height=2, bd=3, relief=tk.SUNKEN)
    self.mwmLabel = tk.Label(self.mouseWheelFrame, text="Mouse wheel manipulation")
    self.mwmLabel.pack()
    # Modify by multiplication or addition (modModus)
    self.modModus = tk.StringVar()
    self.modModus.set("mul")
    # Entry field and radiobutton to specify factor to be used
    self.factorFrame = tk.Frame(self.mouseWheelFrame)
    self.modEntryTextMul = tk.StringVar()
    self.modEntryFactor = tk.Entry(self.factorFrame, textvariable=self.modEntryTextMul, width=6)
    self.modEntryFactor.pack(side=tk.LEFT)
    self.radioMultipli = tk.Radiobutton(self.factorFrame, text="Multiply", value="mul", variable=self.modModus)
    self.radioMultipli.pack(side=tk.LEFT)
    self.factorFrame.pack(fill=tk.BOTH)
    # Entry field and radiobutton to specify step (delta) to be used
    self.addFrame = tk.Frame(self.mouseWheelFrame)
    self.modEntryTextAdd = tk.StringVar()
    self.modEntryAdd = tk.Entry(self.addFrame, textvariable=self.modEntryTextAdd, width=6)
    self.modEntryAdd.pack(side=tk.LEFT)
    self.radioAdd = tk.Radiobutton(self.addFrame, text="Add", value="add", variable=self.modModus)
    self.radioAdd.pack(side=tk.LEFT)
    self.addFrame.pack(fill=tk.BOTH)
    # Set text fields for modification factor/step to default
    self.modEntryTextMul.set(self.modProps[self.selectedPar.get()]["modValMul"])
    self.modEntryTextAdd.set(self.modProps[self.selectedPar.get()]["modValAdd"])
    self.modEntryTextAdd.trace("w", self._modModeChangedAdd)
    self.modEntryTextMul.trace("w", self._modModeChangedMul)
    # Show the frame
#     self.mouseWheelFrame.grid(row=4, column=0, columnspan=3, pady=10)
    self.mouseWheelFrame.pack()
    
    
    # React to change in modify Modus
    self.modModus.trace('w', self._modModusChanged)
    # React to a change in the active parameter
    self.selectedPar.trace("w", self._activeParameterChanged)
    
    dummyLabel = tk.Label(self.controlFrame)
    dummyLabel.pack()
    
    # Fit button and fit ranges
    self.fitRangeFrame = tk.Frame(self.controlFrame, bd=3, relief=tk.SUNKEN)
    self.fit_lo = tk.DoubleVar()
    self.fit_hi = tk.DoubleVar()
    self.fit_lo.set(min(plotter.x))
    self.fit_hi.set(max(plotter.x))
    self.fitRangeLoLim = tk.Entry(self.fitRangeFrame, textvariable=self.fit_lo, width=9, bd=2)
    self.fitRangeHiLim = tk.Entry(self.fitRangeFrame, textvariable=self.fit_hi, width=9, bd=2)
    self.fitRangeLabel = tk.Label(self.fitRangeFrame, text="Fit range:")
    self.fitButton = tk.Button(self.fitRangeFrame, text="Fit", command=self._fitClicked)
    self.fitButton.pack(side=tk.BOTTOM, fill=tk.X)
    self.fitRangeLabel.pack(side=tk.LEFT, fill=tk.X)
    self.fitRangeHiLim.pack(side=tk.RIGHT)    
    self.fitRangeLoLim.pack(side=tk.RIGHT)
    self.fitRangeFrame.pack(fill=tk.X)
    #self.fitRangeLoLim.bind('<Return>', self._fitRangeChanged())
    self.fit_lo.trace("w", self._fitRangeChanged)
    self.fit_hi.trace("w", self._fitRangeChanged)
    self.numberClicked=0
    #self.modModus.trace('w', self._modModusChanged)
    #self.modModus.trace('w', self._modModusChanged)
    
    dummyLabel = tk.Label(self.controlFrame)
    dummyLabel.pack()
    
    self.showFrame = tk.Frame(self.controlFrame, bd=3)#, relief=tk.SUNKEN)
    self.showFrame.pack(side=tk.TOP)
    
    self.parSumButton = tk.Button(self.showFrame, text="Parameter summary", command=self._parameterSummaryClicked)
    self.parSumButton.pack(side=tk.LEFT)
        
    self.valSetButton = tk.Button(self.showFrame, text="Value set code", command=self._valueSetClicked)
#     self.valSetButton.grid(row=7, column=2)
    self.valSetButton.pack(side=tk.RIGHT)
     
    # a tk.DrawingArea
#     self.canvas.get_tk_widget().grid(column=0, columnspan=7, row=0, rowspan=10)
    self.canvas.get_tk_widget().pack()
    self.cid = self.f.canvas.mpl_connect('button_press_event', self._mouseButtonClicked)
    self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.plotFrame)
    self.toolbar.update()
    self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _quit():
      # stops main loop
      self.root.quit()
      # this is necessary on Windows to prevent
      # Fatal Python Error: PyEval_RestoreThread: NULL tstate
      self.root.destroy()

    self.quitButton = tk.Button(master=self.controlFrame, text='Quit', command=_quit)
    self.quitButton.pack(side=tk.BOTTOM)
    
    # Plot the model for the first time
    self._parameterValueChanged()
    
    # Whether or not parameter summary exists
    self.root.parSumWin = None
    
    if platform.system() == "Darwin":
        self.root.bind("<MouseWheel>", self._onWheel) # for OS X

   
  def _parameterActivated(self, *args):
      print(args)

  def _mouseButtonClicked(self, event):
    #print "mouse Button Clicked: ",event
    if event.button != 2: return
    # Check whether click occurred on plot axes.
    if event.inaxes is not self.plotter.a: 
      #print "not in axis"
      #print event.inaxes, self.plotter.a
      return
    else: 
      #print event.xdata, event.ydata
      xrng = self.plotter.a.get_xlim()
      #print 
      digits = int(numpy.log10(1./((xrng[1] - xrng[0])/1000.))+2)
      #print digits
      if self.numberClicked == 0:
          self.fit_lo.set(round(event.xdata,digits))
          self.numberClicked=1
      else:
          self.fit_hi.set(round(event.xdata, digits))
          self.numberClicked=0
      self.plotter.plot(self.f, self.odf)
      
      self._fitRangeChanged()

  def _fitRangeChanged(self, *args):
      #print "fitRangeChanged", args
      #print self.fit_lo.get(), self.fit_hi.get()
      
      if self.leftMask:
          self.plotter.a.collections.remove(self.leftMask)
      if self.rightMask:    
          self.plotter.a.collections.remove(self.rightMask)
      if self.activeLine:
          self.activeLine[0].remove()
      try:
        y0, y1 = self.plotter.a.get_ylim()
        y0, y1 = y0*10, y1*10
        if y0 > 0: y0 = 0
        if y1 < 0: y1 = 0.
        self.leftMask = self.plotter.a.fill_between([min(self.plotter.x),self.fit_lo.get()],[y0,y0],y2=[y1,y1],alpha=0.2)
        self.rightMask = self.plotter.a.fill_between([self.fit_hi.get(), max(self.plotter.x)],[y0,y0],y2=[y1,y1],alpha=0.2)
      except:
          pass
      #self._parameterValueChanged()
      self.plotter.chi2 = self.plotter.get_chi2(self.odf, lims=[self.fit_lo.get(), self.fit_hi.get()])
      self._updateChi2()
      self._updateDof()
      #print self.chi2value.get(), self.fit_lo.get(), self.fit_hi.get()
      self.f.canvas.draw()    

  def _updateDof(self):
      self.dofValue.set(len(self.plotter.fitIdx)- len(self.odf.freeParameters())-1)
      
  def _updateChi2(self):
      try:
        if self.plotter.chi2 > 1e-3 and self.plotter.chi2 < 1e4:
          self.chi2value.set("%8.3f" % self.plotter.chi2)
        else:
          self.chi2value.set("%8.3e" % self.plotter.chi2)
      except:
        self.chi2value.set("N/A")
  
  def _parameterSummaryClicked(self, *args):
    """
    """
    if self.root.parSumWin is None:
      self.root.parSumWin = ShowParameterSummary(self.root, self.odf)
    self.root.parSumWin.updateInfo(self.odf)

  def _valueSetClicked(self, *args):
    """
    """
    if self.root.parSumWin is None:
      self.root.parSumWin = ValueSetSummary(self.root, self.odf)
    self.root.parSumWin.updateInfo(self.odf)

  def _modModusChanged(self, *args):
    """
      Modus for modification changed (multiply or add)
    """
    self.modProps[self.selectedPar.get()]["modus"] = self.modModus.get()

  def _thawFrozenChange(self, *args):
    """
      Called on click to freeze/thaw radiobuttons
    """
    if self.parIsFree.get() == True:
      self.odf.thaw(self.selectedPar.get())
    else:
      self.odf.freeze(self.selectedPar.get())
    self._updateDof()

  def _setToClicked(self):
    """
      Called when "Set to" button is hit
    """
    x = SetToDialog(self.root, self.odf[self.selectedPar.get()], self.selectedPar.get())
    if x.newVal is not None:
      self.odf[self.selectedPar.get()] = x.newVal
      self._parameterValueChanged()

  def _fitClicked(self):
    """
    """
    self.plotter.fit(self.odf, lims=[self.fit_lo.get(), self.fit_hi.get()])
    self._parameterValueChanged(allFree = True)

  def _modModeChangedAdd(self, *args):
      self.modProps[self.selectedPar.get()]["modValAdd"] = self.modEntryTextAdd.get()
      self.modProps[self.selectedPar.get()]["modus"] = "add"
      self.modModus.set("add")
      
  def _modModeChangedMul(self, *args):
      self.modProps[self.selectedPar.get()]["modValMul"] = self.modEntryTextMul.get()
      self.modProps[self.selectedPar.get()]["modus"] = "mul"
      self.modModus.set("mul")

  def _modModeChanged(self, *args):
    """
      Called when mode of modification (add/mul)  is changed
    """
    self.modProps[self.selectedPar.get()]["modValAdd"] = self.modEntryTextAdd.get()
    self.modProps[self.selectedPar.get()]["modValMul"] = self.modEntryTextMul.get()

  def _activeParameterChanged(self, *args):
    """
      Called when the currently active parameter (not its value) is changed.
    """
    # Set the control panels back to that parameter's settings
    pname = self.selectedPar.get()
    newText = "Value:  % g" % (self.odf[pname])
    self.modModus.set(self.modProps[pname]["modus"])
    # Take care of multiply/add radiobuttons
    tmp = self.modProps[pname]["modus"]
    self.modEntryTextAdd.set(self.modProps[pname]["modValAdd"])
    self.modEntryTextMul.set(self.modProps[pname]["modValMul"])
    self.modModus.set(tmp)
    if self.modModus.get() == "mul":
      self.radioMultipli.select()
    else:
      self.radioAdd.select()
    #print  
    
  def _parameterValueChanged(self, allFree=False):
    """
      Called when the value of the current parameter is changed.
    """
    # Update value in label and plot new model
    if not allFree:
      newText = "% g" % (self.odf[self.selectedPar.get()])
      #self.valLabel.config(text=newText)
      #print self.selectedPar.get()
      #print "selected Par: ",self.selectedPar.get()," value: ",self.odf[self.selectedPar.get()]
      self.singleParameterVar[self.selectedPar.get()].set(newText)
    else:
      for p in self.odf.freeParamNames():
        newText = "% g" % (self.odf[p])
        self.singleParameterVar[p].set(newText)

    self.plotter.plot(self.f, self.odf)
    #print self.plotter.chi2
    self._updateChi2()
    self.f.canvas.draw()

  def _onWheel(self, event):
    val = self.odf[self.selectedPar.get()]
    pname = self.selectedPar.get()
    try:
      if self.modModus.get() == "add":
        mf = float(self.modEntryTextAdd.get())
      elif self.modModus.get() == "mul":
        mf = float(self.modEntryTextMul.get())
    except ValueError:
      tkMessageBox.showwarning("Invalid float", "Cannot convert " + self.modEntry.get() + " to float." + \
                               " Make it a valid value to proceed.")
      return
    if self.modModus.get() == "mul":
        if event.delta>0:
          self.odf[pname] = val * mf
        else:
          self.odf[pname] = val / mf
    else:
        self.odf[pname] = val + mf * event.delta/2.
    self._parameterValueChanged()

  def _mouseWheel(self, event):
    """
      Mouse wheel moved
    """
    # event.num == 4 -> up
    # event.num == 5 -> down
    #print "mouse wheel  ",event
    val = self.odf[self.selectedPar.get()]
    pname = self.selectedPar.get()
    try:
      if self.modModus.get() == "add":
        mf = float(self.modEntryTextAdd.get())
      elif self.modModus.get() == "mul":
        mf = float(self.modEntryTextMul.get())
    except ValueError:
      tkMessageBox.showwarning("Invalid float", "Cannot convert " + self.modEntry.get() + " to float." + \
                               " Make it a valid value to proceed.")
      return
    if event.num == 4:
      if self.modModus.get() == "mul":
        self.odf[pname] = val * mf
      else:
        self.odf[pname] = val + mf
    elif event.num == 5:
      if self.modModus.get() == "mul":
        self.odf[pname] = val / mf
      else:
        self.odf[pname] = val - mf
    self._parameterValueChanged()

  def show(self):
    """
      Show the GUI.
    """
    self.canvas.show()
    tk.mainloop()


    
# Set the standard for FFModelExplorer    
FFModelExplorer = FFModelExplorerList

def ffmodelExplorer(odf, plotter, version="list"):
  """
    Instantiate the model explorer.
    
    Parameters
    ----------
    odf : Instance of OneDFit
        The model to be adapted.
    plotter : Instance of FFModelPlotFit or custom
        Class instance managing the actual plotting.
    version : string, {list}
        The version of model explorer. Currently, only
        'list' is supported.
    
    Returns
    -------
    ffmod : Model explorer
        An instance of the model explorer.
  """
  if version == "list":
    ffmod = FFModelExplorerList(odf, plotter)
    return ffmod
  elif version == "dropdown":
    PE.warn(PE.PyADeprecationError("Please note that the dropdown version is no longer supported. " + \
                                   "The list version will be called instead."))
    ffmod = FFModelExplorerList(odf, plotter)
    return ffmod
  else:
    raise(PE.PyAValError("Unknown version: '" + str(version) + "'", \
                         where="ffmodelExplorer", \
                         solution="Choose between 'list' and 'dropdown'."))