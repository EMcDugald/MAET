#!/opt/anaconda3/bin/python
#---------------------------------------------------
#
#    dview.py
#
#   A program for displaying and comparing rectangular arrays as a colormap
#   The arrays are read from  .d files
#
#---------------------------------------------------------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
import os
import sys
import utils
# import winsound
#import msvcrt
import time
from tkinter import simpledialog
from tkinter import messagebox

from tkinter import filedialog
from matplotlib.colors import ListedColormap

#------------------------ Subroutines --------------------------
def get_maxplot():
    global  custommax, dmaxplot, data_to_use

    if (custommax == False):
       dmax = np.amax(data_to_use)
       dmin = np.amin(data_to_use)

       if(dmax <= 0.0):
          dmaxplot = - dmin
       else:
          dmaxplot = dmax
       if(dmaxplot < 1.e-30):
          dmaxplot = 1.


def draw_image():

    global data,fig,ax,displayname,cb, sign, dmaxplot,data1, nfiles
    global data_to_use, mymap, title_color, twosided, need_new_figure, data2, datalinf, localpath

    startdraw = time.time()

    title = displayname


    dmax = np.amax(data_to_use)
    dmin = np.amin(data_to_use)

    #print(dmaxplot)

    norm = mpl.colors.Normalize(vmin=-dmaxplot*twosided,vmax=dmaxplot)
    #norm = mpl.colors.Normalize(vmin= 0,vmax=dmaxplot)

    #print(dmax)

    #mpl.rcParams['axes.edgecolor'] = 'black'






    #pos = ax.imshow(data_to_use, origin='lower',cmap=mymap, norm = norm)
    #cb = fig.colorbar(pos, ax = ax)

    try:
       cb.remove()
    except:
       pass

    plt.clf()


    if(need_new_figure):
       need_new_figure = False
       plt.close()

       #print('open new figure')
       fig = plt.figure(1,figsize=(15,9))
       fig.canvas.mpl_connect('key_press_event', on_press)
       plt.show()

    plt.imshow(data_to_use, origin='lower',cmap=mymap, norm = norm)
    cb = plt.colorbar()

    ny, nx = np.shape(data_to_use)

    #firstline  = bytes.decode(str(ny).encode()) + ' x ' + bytes.decode(str(nx).encode()) + '     '+title + '\n \n'

    firstline  = '{:6d}'.format(ny) + ' x ' + '{:6d}'.format(nx) + '           '+title + '\n \n'

    thirdline =  'Max = '  + '{:.4g}'.format(dmax) + '      Min = '  + '{:.4g}'.format(dmin)
    #thirdline =  'Range =   ['  + '{:.4g}'.format(dmax) + ','+'{:.4g}'.format(dmin)+']'

    if(nfiles != 2):
       longtitle = firstline + thirdline + '\n'

    if(nfiles == 2):
       diffl2 = np.sum(data_to_use*data_to_use)
       errl2 = np.sqrt(diffl2/datal2)
       difflinf = np.max(np.abs(data_to_use))
       errlinf = difflinf/datalinf

       longtitle = firstline + thirdline + '         $\Delta_{{L_2}}$  = ' + '{:.3g}'.format(100.*errl2)  + '%          $\Delta_\infty$ = ' + '{:.3g}'.format(100.*errlinf) + '% \n'

    if(sign == -1):
       longtitle = '(-)    ' + longtitle

    plt.title(longtitle, color=title_color)
    #plt.title(longtitle, color="brown",size = 'x-large')

    fig.canvas.draw()
    fig.canvas.flush_events()
    enddraw = time.time()

    print('draw time = ', enddraw - startdraw)

def next_prev(key):
    global data, ax,cb, sign, custommax, dmaxplot, data1, data_to_use, path, name, extension, displayname, lock

    success = True
    lastfour = name[-4:]
    #print(lastfour)
    if(lastfour[-3] in ["r","R","i","I","a","A"]):
       letter = lastfour[-3]
       have_letter = True
       #print(letter)
       stri = lastfour[0] + lastfour[-2:]
    else:
       have_letter = False
       stri = lastfour

    try:
       current_number = int(stri)
       have_number = True
    except:
       # winsound.Beep(int(2000),200)
       have_number = False
       print('---------------------------------------------------------------------')
       print('Cant figure out the name of the next or previous file in the sequence')
       print('---------------------------------------------------------------------')
       success = False

    if(have_number):

       oldname = name

       if (key == 'n'):
          new_number = current_number + 1
       if (key == 'p'):
          new_number = current_number - 1

       if(new_number < 0):
          # winsound.Beep(int(2000),200)
          success = False
       else:
          new_stri = str(new_number).zfill(4)
          #print(new_stri)
          #print(have_letter)
          if(have_letter):
             newname = name[:-4] + new_stri[1]+letter+new_stri[-2:]
          else:
             newname = name[:-4] + new_stri

          newpath = path+newname + extension

          #print(newpath)

          #if(os.path.isfile(newpath)):
          try:
             data_to_use = utils.readd(newpath)
             if(nfiles == 2):
                data_to_use = data_to_use - data1

             if(sign == -1):
                data_to_use = - data_to_use

             name = newname
             (displayname,ext) = os.path.splitext(newpath)

             if(nfiles == 2):
                displayname = displayname + " - " + path_and_name1

             #print(displayname)

             get_maxplot()

             draw_image()
          except:
             # winsound.Beep(int(2000),200)
             success = False
    return success

def on_press(event):
 global data, ax,cb, sign, custommax, dmaxplot, data1, data_to_use, path, name, extension, displayname, lock, style_number, mycolors,mymap,title_color,twosided, need_new_figure

 if (event.key in ['escape','q','Q']):
    quit()

 if(not lock):
    lock = True
    #print('press', event.key)
    fig.canvas.flush_events()


    if event.key == ' ':
       data_to_use = - data_to_use
       sign = - sign

       draw_image()

    if ( event.key in ['m','M']):
       #print('Swithcing mapping')

       if(custommax == True):
          custommax = False
          get_maxplot()
          draw_image()
       else:
          newvalue = simpledialog.askfloat("Rescale", "Please, please, enter new maximum#", initialvalue=dmaxplot)
          if(newvalue != None):
             custommax = True
             dmaxplot = newvalue
             draw_image()

    if ( (event.key == 'n') | (event.key == 'p') ):
       success = next_prev(event.key)
       # winsound.Beep(int(300),80)

    if ( event.key == 'N'):
       success = True
       while (success):
          success = next_prev('n')
          #winsound.Beep(int(300),80)

    if ( event.key == 'P'):
       success = True
       while (success):
          success = next_prev('p')
          #winsound.Beep(int(300),80)

    if( (event.key == 'h') | (event.key == 'H')):
        stri = messagebox.showinfo("Help",helpstring)

    if ( event.key in ['c','b','C','B']):
       style_number = style_number + 1
       if(style_number > 2):
          style_number = 0
       #print("style = ",style_number)
       need_new_figure = True

       plt.style.use("dark_background")

       if(style_number == 0):
          params = {"text.color" : "grey",
          "xtick.color" : "grey",
          "ytick.color" : "grey",
          "axes.edgecolor" : "black",
          "patch.edgecolor" : "grey",
          "savefig.directory" : localpath,
          "savefig.format" : "pdf"}


          plt.rcParams.update(params)

          mymap = ListedColormap(mycolors)
          title_color = "yellow"

       if(style_number > 0):
          plt.style.use("default")

          params = {"xtick.color" : "black",
          "ytick.color" : "black",
          "axes.edgecolor" : "white",
          "savefig.directory" : localpath,
          "savefig.format": "pdf"}



          plt.rcParams.update(params)

          mymap = plt.cm.get_cmap('gist_gray')
          title_color = "black"

       if(style_number == 2):
          twosided = 0
       else:
          twosided = 1

       draw_image()

 lock = False

#------------------- Start main program -----------------------------



arguments = len(sys.argv) - 1


print("Starting dview")
print("**************\n")

#print('Argument 1 ',sys.argv[1])

#time.sleep(10)


root = tk.Tk()
#root.geometry("0x0+0+0")
root.withdraw()

if(arguments not in [0, 1, 2]):
    print("nunber of arguments is ",arguments," WTF ???")
    exit()

if(arguments == 0):
   # Dialog to select .d file
   root.withdraw()
   fullpath = filedialog.askopenfilename(filetypes=[("Data files", "*.d*")])

   if not fullpath:
      exit()
   nfiles = 1

elif(arguments == 1):
   fullpath=sys.argv[1]
   nfiles = 1

else:
   fullpath =sys.argv[1]
   fullpath1=sys.argv[2]

(path_and_name,extension) = os.path.splitext(fullpath)
(path,name) = os.path.split(path_and_name)

if(len(path) > 0):
   path = path+'\\'
#print(path)
#print(name)
#print(extension)

if(arguments ==2 ):
   nfiles = 2
   (path_and_name1,extension1) = os.path.splitext(fullpath1)
   (path1,name1) = os.path.split(path_and_name1)


   if(len(path1) > 0):
      path1 = path1+'\\'

   #print(path1)
   #print(name1)
   #print(extension1)
#----------------------------------- now let's start reading ----------------------------

data = utils.readd(fullpath)

data_to_use = data

displayname = path_and_name

if(nfiles == 2):
   data1 = utils.readd(fullpath1)

   if(np.shape(data1) != np.shape(data)):
       print('Data dimensions incompatible')
       quit()

   data_to_use = data - data1    # this needs to be changed

   displayname = displayname + " - " + path_and_name1


datal2   = np.sum(data*data)
datalinf = np.max(np.abs(data))


#------------- set global parameters ---------------------------
#print(os.path.splitext(fpath))
localpath = os.getcwd()
#print(localpath)

sign = 1
lock = False

mpl.rcParams["savefig.directory"] = localpath
mpl.rcParams["savefig.format"] = "pdf"


plt.style.use("dark_background")

fig = plt.figure(1,figsize=(15,9))
fig.canvas.mpl_connect('key_press_event', on_press)

params = {"text.color" : "grey",
"xtick.color" : "grey",
"ytick.color" : "grey",
"axes.edgecolor" : "black",
"patch.edgecolor" : "grey",
"figure.edgecolor" : "red"}


plt.rcParams.update(params)

style_number = 0
title_color = "yellow"
twosided = 1
need_new_figure = False

mycolors=utils.mycolormap()
mymap = ListedColormap(mycolors)


custommax = False

helpstring = "List of the keybord shortcuts:  \n"                              \
+"==============================  \n"                                          \
+"      \n"                                                                \
+" h   =   Show help                                     \n"\
+" n   =   Show the next file in a numbered sequence     \n"\
+" P   =   Show the previous file in a numbered sequence \n"\
+" N   =   Play a numbered sequence of fiules forward    \n"\
+" N   =   Play a numbered sequence of fiules backwards  \n"\
+" p   =   Show the previous file in a numbered sequence \n"\
+" c|b =   Cycle through available color schemes         \n"\
+" m   =   Manual scaling of the colors                  \n"\
+"' '  =   Switching the sign of the data                  \n"\
+" s   =   Save image as a .pdf                          \n"\
+" q   =   Quit                                          \n"\
+"     "

print(helpstring)



get_maxplot()

draw_image()

plt.show()

quit()






#plt.rcParams.update({
#    "lines.color": "white",
#    "patch.edgecolor": "white",
#    "text.color": "black",
#    "axes.facecolor": "white",
#    "axes.edgecolor": "lightgray",
#    "axes.labelcolor": "white",
#    "xtick.color": "white",
#    "ytick.color": "white",
#    "grid.color": "lightgray",
#    "figure.facecolor": "black",
#    "figure.edgecolor": "black",
#    "savefig.facecolor": "black",
#    "savefig.edgecolor": "black"})
