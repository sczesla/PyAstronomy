import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from PyAstronomy.pyaC import pyaErrors as PE

import Tkinter as tk


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
    """
    return "%g, %g" % (self.xdata, self.ydata)


class Picker:
  """
    GUI for interactive point selection.
    
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
  
  def __init__(self):
    self.windowTitle = "PyA Picker"
    self.f = Figure()
    self.a = self.f.add_subplot(111)
    
    # Active and Inactive plot style
    self.astyle = "ro"
    self.istyle = "yp"
    
    self.pointList = []

    self.root = tk.Tk()
    # Make the widgets expand/shrink as window size changes
    self.root.columnconfigure(0, weight=1)
    self.root.rowconfigure(0, weight=1)
    
    # A frame containing the mpl plot
    self.plotFrame = tk.Frame()
    self.plotFrame.grid(column=0, columnspan=7, row=0, rowspan=10, sticky="nsew")
    self.canvas = FigureCanvasTkAgg(self.f, master=self.plotFrame)
    
    # A frame containing the box with selected points
    # and control buttons
    self.pointFrame = tk.Frame(self.root)
    self.pointFrame.grid(column=7, columnspan=3, row=0, rowspan=10)
    
    self.lb = tk.Listbox(self.pointFrame, height=15)
    self.lb.grid(row=0,rowspan=6, column=0, columnspan=1)
    self.lb.bind("<<ListboxSelect>>", self._lbSelect)
    
    self.copyButton = tk.Button(self.pointFrame, text="Open List", command=self._showList)
    self.copyButton.grid(row=8, sticky='w')
    
    self.removeButton = tk.Button(master=self.pointFrame, text="Remove", \
                                  command=self._removeButtonClicked)
    self.removeButton.grid_configure(row=7, sticky="w")

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

    self.quitButton = tk.Button(master=self.pointFrame, text='Quit', command=_quit)
    self.quitButton.grid_configure(row=9, sticky="e", pady=30)
    
  def _showList(self):
      """
        Creates a new window containing a copy and paste aware list with the 
        points.
      """
      def _quitWin():
          win.destroy()
      win = tk.Tk()
      text =  tk.Text(master=win)
      text.pack()
      for p in self.pointList:
        text.insert(tk.END, p.asStr()+"\n")
      quit = tk.Button(master=win, text='Quit', command=_quitWin)
      quit.pack()
      win.wm_title("Copy Window")

    

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
      pli, lli = self._seachPoint(point.lbIdent)
    else:
      pli, lli = self._seachPoint(lbString)
    # Remove 'old' point from plot
    self.a.lines.pop(lli)
    if state == "active":
      style = self.astyle
    elif state == "inactive":
      style = self.istyle
    # Add new point in specified color
    self.a.plot([self.pointList[pli].xdata], [self.pointList[pli].ydata], style)
    self.pointList[pli].mplLine = self.a.lines[-1]
      
  def _seachPoint(self, lbString):
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
      React on change of listbox selection.
      
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
    self.f.canvas.draw()
    
    # Insert the point's content into the list box
    self.lb.insert(tk.END, p.asStr())
    # Save the "list box identifier"
    p.lbIdent = p.asStr()
    # Save the point
    self.pointList.append(p)

  def _removeButtonClicked(self):
    """
      Remove a point.
    """
    # Check whether anything was selected
    sel = self.lb.curselection()
    if len(sel) == 0: return
    # Find indices for the point 
    el = self.lb.get(int(sel[0]))
    pli, lli = self._seachPoint(el)
    # Remove it
    self.a.lines.pop(lli)
    self.pointList.pop(pli)
    self.lb.delete(sel[0])
    self._updateView()

  def _updateView(self):
    """
      Redraw MPL canvas
    """
    self.f.canvas.draw()

  def pick(self):
    """
      Select points.
      
      .. note:: Before selection points, a plot should
                have been created from which points may
                be selected.
      
      Returns
      -------
      Points : List of tuples
          A list of two-element tuples, holding the x and
          y coordinates of the selected points.
    """
    self.root.wm_title(self.windowTitle)
    self.canvas.show()
    tk.mainloop()
    # Prepare return value
    points = []
    for p in self.pointList:
      points.append((p.xdata, p.ydata))
    return points
    
    
