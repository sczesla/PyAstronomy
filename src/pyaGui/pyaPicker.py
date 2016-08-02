import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from PyAstronomy.pyaC import pyaErrors as PE

import six.moves.tkinter as tk
import pickle


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


class Picker:
  """
    GUI for interactive point selection.
    
    Parameters
    ----------
    saveFile : string, optional
        If a filename is given, the list of points will be
        saved to a file with that name using pickle when the
        class instance if destructed. The points are saved as
        a python list of two-float-tuples, holding the x and y
        coordinates of the points. 
    
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
  
  def __del__(self):
    """
      Destruct instance
    """
    if not self.saveFile is None:
      # Save current point list to file if filename was given
      pickle.dump(self.points, open(self.saveFile, 'w'))
  
  def __init__(self, saveFile=None):
    self.windowTitle = "PyA Picker"
    self.f = Figure()
    self.a = self.f.add_subplot(111)
    
    # Store save filename used to dump data on
    # destruction of class instance
    self.saveFile = saveFile
    
    # Active and Inactive plot style
    self.astyle = "ro"
    self.istyle = "yp"
    
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
    
    self.lb = tk.Listbox(self.pointFrame, height=15)
    self.lb.pack(side=tk.TOP)
    self.lb.bind("<<ListboxSelect>>", self._lbSelect)
    
    self.copyButton = tk.Button(self.pointFrame, text="Open List", command=self._showList)
    self.copyButton.pack(side=tk.TOP, fill=tk.X)
    
    self.removeButton = tk.Button(master=self.pointFrame, text="Remove", \
                                  command=self._removeButtonClicked)
    self.removeButton.pack(side=tk.TOP, fill=tk.X)

    # a tk.DrawingArea
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

    self.quitButton = tk.Button(master=self.pointFrame, text='Quit', command=_quit)
    self.quitButton.pack(side=tk.BOTTOM, fill=tk.X)
    
  def _showList(self):
      """
        Creates a new window containing a copy-and-paste aware list of the points.
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
    pli, lli = self._searchPoint(el)
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
    self.points = []
    for p in self.pointList:
      self.points.append((p.xdata, p.ydata))
    return self.points[:]
    
    
