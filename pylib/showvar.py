import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot
from matplotlib import cm
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.figure import Figure


import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

root = Tk.Tk()
root.wm_title("Embedding in TK")


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate


def on_key_event(event):
    print('you pressed %s' % event.key)
    key_press_handler(event, canvas, toolbar)




def drawColorBar(data, CS):

    global colorBarPresent, cbar, a, f, cax

    cbar = f.colorbar(CS, cax=cax)

    #-number of ticks on colorbar
    nlevel=5

    vmin=np.min(data)
    vmax=np.max(data)
    cbar_ticks = np.linspace(vmin, vmax, num=nlevel, endpoint=True)
    cbar.set_ticks(cbar_ticks)

    cbar.update_normal(CS)
    cbar.draw_all()
    cbar.update_bruteforce(CS)





def MakePlot(data, xdata, ydata, title, xtitle, ytitle, showcanvas=True,):

    global colorBarPresent, cbar, a, f, cax

    cax.clear()
    a.clear()

    ncolor=20  #150 #-number of color levels
    CS = a.contourf(xdata, ydata, data, ncolor, cmap=cm.jet)
    a.set_title(title)
    a.set_xlabel(xtitle)
    a.set_ylabel(ytitle)
    a.grid(True)
    
    drawColorBar(data, CS)
    if showcanvas:
        canvas.show()




def _showte(showcanvas=True):
    ##print('entering showte()')
    ev=1.6022e-19
    data = bbb.te[:,:]/ev
    xdata=com.rm[:,:,0]
    ydata=com.zm[:,:,0]
    title='T$\mathregular{_e}$ [eV]'
    xtitle='R [m]'
    ytitle='Z [m]'
    MakePlot(data, xdata, ydata, title, xtitle, ytitle, showcanvas=showcanvas)
    ##print('leaving showte()')



def _showti():
    ##print('entering showti()')
    ev=1.6022e-19
    data = bbb.ti[:,:]/ev
    xdata=com.rm[:,:,0]
    ydata=com.zm[:,:,0]
    title='T$\mathregular{_i}$ [eV]'
    xtitle='R [m]'
    ytitle='Z [m]'
    MakePlot(data, xdata, ydata, title, xtitle, ytitle)
    ##print('leaving showti()')



def _showni():
    ##print('entering showni()')
    data = bbb.ni[:,:,0]
    xdata=com.rm[:,:,0]
    ydata=com.zm[:,:,0]
    title='N$\mathregular{_i}$ [m-3]'
    xtitle='R [m]'
    ytitle='Z [m]'
    MakePlot(data, xdata, ydata, title, xtitle, ytitle)
    ##print('leaving showni()')


def _showng():
    ##print('entering showni()')
    data = bbb.ng[:,:,0]
    xdata=com.rm[:,:,0]
    ydata=com.zm[:,:,0]
    title='N$\mathregular{_g}$ [m-3]'
    xtitle='R [m]'
    ytitle='Z [m]'
    MakePlot(data, xdata, ydata, title, xtitle, ytitle)
    ##print('leaving showni()')


def _showup():
    ##print('entering showni()')
    data = bbb.up[:,:,0]
    xdata=com.rm[:,:,0]
    ydata=com.zm[:,:,0]
    title='U$\mathregular{_{||i}}$ [m/s]'
    xtitle='R [m]'
    ytitle='Z [m]'
    MakePlot(data, xdata, ydata, title, xtitle, ytitle)
    ##print('leaving showni()')







########################################################################

f = Figure(figsize=(10, 8), dpi=100)
a = f.add_subplot(111)
cax = make_axes_locatable(a).append_axes("right", size="5%", pad=0.05)
colorBarPresent=False

#-make initial plot of Te[eV]
_showte(showcanvas=False)


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

canvas.mpl_connect('key_press_event', on_key_event)



button1 = Tk.Button(master=root, text='Show Te', command=_showte)
button1.pack(side=Tk.LEFT)

button2 = Tk.Button(master=root, text='Show Ti', command=_showti)
button2.pack(side=Tk.LEFT)

button3 = Tk.Button(master=root, text='Show Ni', command=_showni)
button3.pack(side=Tk.LEFT)

button4 = Tk.Button(master=root, text='Show Up', command=_showup)
button4.pack(side=Tk.LEFT)

button5 = Tk.Button(master=root, text='Show Ng', command=_showng)
button5.pack(side=Tk.LEFT)

button = Tk.Button(master=root, text='Quit', command=_quit)
button.pack(side=Tk.LEFT)



Tk.mainloop()
# If you put root.destroy() here, it will cause an error if
# the window is closed with the window manager.
