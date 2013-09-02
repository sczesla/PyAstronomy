# -*- coding: utf-8 -*-
import matplotlib

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from PyAstronomy.pyaC import pyaErrors as PE

import Tkinter as tk
import tkMessageBox


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
      self.a.lines.pop(-1)
    # Plot the model
    self.modelLine = self.a.plot(self.x, model, 'r--')
    
    if self.withResiduals:
      if self.residualLine is not None:
        self.ar.lines.pop(-1)
      res = self.y - model
      if self.yerr is not None:
        res /= self.yerr
      self.residualLine = self.ar.plot(self.x, res, 'bp')
    if firstCall:
      self.a.autoscale(False)

  def fit(self, odf):
    """
      Fit the model.
      
      Parameters
      ----------
      odf : funcFit model class
          The model to be fitted.
    """
    odf.fit(self.x, self.y, yerr=self.yerr)


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
    self.showText.config(width=max(map(lambda x: len(x), lines)), \
                         height=len(lines)+1)
    t = "".join(map(lambda x: x + "\n", lines))
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
    
    self.showText.config(width=max(map(lambda x: len(x), lines)), \
                         height=len(lines)+1)
    t = "".join(map(lambda x: x + "\n", lines))
    self.showText.insert(tk.END, t)
    self.showText.config(state=tk.DISABLED)

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
  """
  
  def __init__(self, odf, plotter):
    self.f = Figure()
    
    # Save reference to the model
    self.odf = odf
    # Save reference to the plotter
    self.plotter = plotter
    
    self.root = tk.Tk()
    self.root.wm_title("PyA Model Explorer")
    # Make the widgets expand/shrink as window size changes
    self.root.columnconfigure(0, weight=1)
    self.root.rowconfigure(0, weight=1)
    
    # Bind the mouse wheel
    self.root.bind("<Button-4>", self._mouseWheel)
    self.root.bind("<Button-5>", self._mouseWheel)
        
    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.grid(column=0, columnspan=7, row=0, rowspan=10, sticky="nsew")
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.controlFrame = tk.Frame(self.root)
    self.controlFrame.grid(column=7, columnspan=3, row=0, rowspan=10, sticky="nsew")
    self.selectedPar = tk.StringVar(self.controlFrame)

    # Get parameters of model
    ps = self.odf.parameters().keys()
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
        return change
    
    # define what happens when a value is set to the entered parameter (closures)
    def hitReturn(k):
      def change(*args):
        self.selectedPar.set(k)
        self.odf[k] = float(self.singleParameterVar[k].get())
        self._parameterValueChanged()
      return change  
    
    # defines what heppens when a parameter's value is changed, but not set yet, i.e., not current (closures)
    def parameterValueChanged(k):
        def valueChanged(*args):
            pp,ll = str(self.singleParameterVar[k].get()).find("."), len(str(self.singleParameterVar[k].get()))
            if round(self.odf[k],ll-pp-1) != round(self.singleParameterVar[k].get(), ll-pp-1):
                self.singleParameterEntry[k].configure(bg = "red")
            else:
                self.singleParameterEntry[k].configure(bg = "white")
        return valueChanged    
    
    # Create an entry for each parameter
    for k in sorted(self.odf.parameters().keys()):
      x = tk.Frame(self.parameterFrame, height=2, bd=2, relief=tk.SUNKEN,pady=2)
      self.singleParameterFrames[k] = x 
      y0 = tk.Radiobutton(x, text=k,variable=self.selectedPar,value=k, width=3)
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
    self.parameterFrame.grid(row=2, column=0, columnspan=3)
    
    # Set of the menu to select the current parameter
    self.selectedPar.set(ps[0])
    
    # Frame to bundle mouse-wheel inputs
    self.mouseWheelFrame = tk.Frame(self.controlFrame, height=2, bd=1, relief=tk.SUNKEN)
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
    self.mouseWheelFrame.grid(row=3, column=0, columnspan=3, pady=10)
    
    # React to change in modify Modus
    self.modModus.trace('w', self._modModusChanged)
    # React to a change in the active parameter
    self.selectedPar.trace("w", self._activeParameterChanged)
    
    # Show the fit button
    self.fitButton = tk.Button(self.controlFrame, text="Fit", command=self._fitClicked)
    self.fitButton.grid(row=5, column=0, columnspan=3)
    
    self.parSumButton = tk.Button(self.controlFrame, text="Parameter summary", command=self._parameterSummaryClicked)
    self.parSumButton.grid(row=6)
        
    self.valSetButton = tk.Button(self.controlFrame, text="Value set code", command=self._valueSetClicked)
    self.valSetButton.grid(row=6, column=2)
     
    # a tk.DrawingArea
    self.canvas.get_tk_widget().grid(column=0, columnspan=7, row=0, rowspan=10)
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
    self.quitButton.grid_configure(row=9, column=2, sticky="se")
    
    # Plot the model for the first time
    self._parameterValueChanged()
    
    # Whether or not parameter summary exists
    self.root.parSumWin = None

  def _parameterActivated(self, *args):
      print args

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
    self.plotter.fit(self.odf)
    self._parameterValueChanged()

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
    print  
    
  def _parameterValueChanged(self):
    """
      Called when the value of the current parameter is changed.
    """
    # Update value in label and plot new model
    newText = "% g" % (self.odf[self.selectedPar.get()])
    #self.valLabel.config(text=newText)
    #print self.selectedPar.get()
    #print "selected Par: ",self.selectedPar.get()," value: ",self.odf[self.selectedPar.get()]
    self.singleParameterVar[self.selectedPar.get()].set(newText)

    self.plotter.plot(self.f, self.odf)
    self.f.canvas.draw()

  def _mouseWheel(self, event):
    """
      Mouse wheel moved
    """
    # event.num == 4 -> up
    # event.num == 5 -> down
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

  def _mouseButtonClicked(self, even):
    pass

  def show(self):
    """
      Show the GUI.
    """
    self.canvas.show()
    tk.mainloop()


class FFModelExplorerDropDownMenu:
  """
    Adapt model interactively.
  """
  
  def __init__(self, odf, plotter):
    self.f = Figure()
    
    # Save reference to the model
    self.odf = odf
    # Save reference to the plotter
    self.plotter = plotter
    
    self.root = tk.Tk()
    self.root.wm_title("PyA Model Explorer")
    # Make the widgets expand/shrink as window size changes
    self.root.columnconfigure(0, weight=1)
    self.root.rowconfigure(0, weight=1)
    
    # Bind the mouse wheel
    self.root.bind("<Button-4>", self._mouseWheel)
    self.root.bind("<Button-5>", self._mouseWheel)
        
    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.grid(column=0, columnspan=7, row=0, rowspan=10, sticky="nsew")
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.controlFrame = tk.Frame(self.root)
    self.controlFrame.grid(column=7, columnspan=3, row=0, rowspan=10, sticky="nsew")
    
    # Get parameters of model
    ps = self.odf.parameters().keys()
    # Set default modification properties
    # Saves these properties for all parameters
    self.modProps = {}
    for p in ps:
      self.modProps[p] = {"modus":"mul", "modValMul":1.02, "modValAdd":0.01}
    
    
    # Bundle "value label" and "set to value" button
    self.valueFrame = tk.Frame(self.controlFrame, height=2, bd=1, relief=tk.SUNKEN)
    # A label to display the variable's value
    self.valLabel = tk.Label(self.valueFrame)
    # Show "Set to" button
    self.setToButton = tk.Button(self.valueFrame, text="Set to value", command=self._setToClicked)
    self.valLabel.pack()
    self.setToButton.pack()
    self.valueFrame.grid(row=1, column=0, columnspan=3)
    
    # Set of the menu to select the current parameter
    self.selectedPar = tk.StringVar(self.controlFrame)
    self.selectedPar.set(ps[0])
    self.pselect = tk.OptionMenu(self.controlFrame, self.selectedPar, *ps)
    self.pselect.grid(row=0,rowspan=1, column=0, columnspan=3, sticky="ew", pady=10)
    
    # Frame to bundle mouse-wheel inputs
    self.mouseWheelFrame = tk.Frame(self.controlFrame, height=2, bd=1, relief=tk.SUNKEN)
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
    self.modEntryTextAdd.trace("w", self._modModeChanged)
    self.modEntryTextMul.trace("w", self._modModeChanged)
    # Show the frame
    self.mouseWheelFrame.grid(row=3, column=0, columnspan=3, pady=10)
    
    # Thaw/Freeze choice
    self.thawFreezeFrame = tk.Frame(self.controlFrame, height=2, bd=1, relief=tk.SUNKEN)
    self.parIsFree = tk.BooleanVar()
    self.radioThaw = tk.Radiobutton(self.thawFreezeFrame, text="Thawed", value=True, variable=self.parIsFree)
    self.radioFreeze = tk.Radiobutton(self.thawFreezeFrame, text="Frozen", value=False, variable=self.parIsFree)
    self.parIsFree.set(self.selectedPar.get() in self.odf.freeParamNames())
    self.radioThaw.pack(side=tk.LEFT)
    self.radioFreeze.pack(side=tk.RIGHT)
    self.thawFreezeFrame.grid(row=4, column=0, columnspan=3, pady=10)
    
    # React to change in modify Modus
    self.modModus.trace('w', self._modModusChanged)
    # React to change in thaw/frozen state
    self.parIsFree.trace('w', self._thawFrozenChange)
    # React to a change in the active parameter
    self.selectedPar.trace("w", self._activeParameterChanged)
    
    # Show the fit button
    self.fitButton = tk.Button(self.controlFrame, text="Fit", command=self._fitClicked)
    self.fitButton.grid(row=5, column=0, columnspan=3)
    
    self.parSumButton = tk.Button(self.controlFrame, text="Parameter summary", command=self._parameterSummaryClicked)
    self.parSumButton.grid(row=6)

    # a tk.DrawingArea
    self.canvas.get_tk_widget().grid(column=0, columnspan=7, row=0, rowspan=10)
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
    self.quitButton.grid_configure(row=9, column=2, sticky="se")
    
    # Plot the model for the first time
    self._parameterValueChanged()
    
    # Whether or not parameter summary exists
    self.root.parSumWin = None

  def _parameterSummaryClicked(self, *args):
    """
    """
    if self.root.parSumWin is None:
      self.root.parSumWin = ShowParameterSummary(self.root, self.odf)
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
    self.plotter.fit(self.odf)
    self._parameterValueChanged()

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
    self.valLabel.config(text=newText)
    self.modModus.set(self.modProps[pname]["modus"])
    # Take care of multiply/add radiobuttons
    if self.modModus.get() == "mul":
      self.radioMultipli.select()
    else:
      self.radioAdd.select()
    self.modEntryTextAdd.set(self.modProps[pname]["modValAdd"])
    self.modEntryTextMul.set(self.modProps[pname]["modValMul"])
    self.parIsFree.set(pname in self.odf.freeParamNames())
    # Take care of frozen/thawed radiobuttons
    if self.parIsFree.get():
      self.radioThaw.select()
    else:
      self.radioFreeze.select()

  def _parameterValueChanged(self):
    """
      Called when the value of the current parameter is changed.
    """
    # Update value in label and plot new model
    newText = "Value:  % g" % (self.odf[self.selectedPar.get()])
    self.valLabel.config(text=newText)
    #print dir(self.plotter.a)
    #print self.plotter.a
    #x, y = self.plotter.a.get_xlim(), self.plotter.a.get_ylim()
    #print x,y
    self.plotter.plot(self.f, self.odf)
    self.f.canvas.draw()

  def _mouseWheel(self, event):
    """
      Mouse wheel moved
    """
    # event.num == 4 -> up
    # event.num == 5 -> down
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

  def _mouseButtonClicked(self, even):
    pass

  def show(self):
    """
      Show the GUI.
    """
    self.canvas.show()
    tk.mainloop()
    
# Set the standard for FFModelExplorer    
FFModelExplorer = FFModelExplorerList