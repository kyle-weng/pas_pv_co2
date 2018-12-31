import os
import pandas
from matplotlib import pyplot
from matplotlib import patches
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy import array

def plotyear(yearselect):
    if yearselect == 2004:
        patchcolor = 'red'
    elif yearselect == 2006:
        patchcolor = 'blue'
    elif yearselect == 2005:
        patchcolor = 'green'
    slopes=[]
    pressure = [1000, 985, 969.9998, 954.9998, 940.0001, 925, 909.9998, 895, 879.9999,      865, 850.0001, 835.0001, 820.0002, 800.0002, 774.9998, 750.0001,      725.0001, 700.0001, 675.0001, 637.5002, 600.0002, 562.5001, 525.0001,      487.5001, 450.0001, 412.5001, 375.0001, 337.5001, 288.0831, 244.875,      208.152, 176.93, 150.393, 127.837, 108.663, 92.36572, 78.51231, 56.38791,      40.17541, 28.36781, 19.7916, 9.292942, 4.076571, 1.65079, 0.6167791,      0.211349, 0.066]
    palt=pandas.read_excel("C:\Users\Kyle W\Desktop\pva.xlsx")
    def slope(x1, y1,x2,y2):
        return (y2 - y1) / (x2 - x1)
    for x in range (4,len(palt["altjan"])):
        slopes.append(slope(1,palt["altjan"][x],7,palt["altjul"][x]))
    domain=range(len(palt["altjan"]))
    domain.remove(0)
    domain.remove(1)
    domain.remove(2)
    domain.remove(3)
    def extrap1d(interpolator):
        xs = interpolator.x
        ys = interpolator.y
    
        def pointwise(x):
            if x < xs[0]:
                return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            elif x > xs[-1]:
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
                return interpolator(x)
    
        def ufunclike(xs):
            return array(map(pointwise, array(xs)))
    
        return ufunclike
    f_i = interp1d(domain, slopes)
    f_x = extrap1d(f_i)
    
    slopesextended=list(f_x([0,1,2,3]))
    for x in range (0,len(slopesextended)):
        palt["altjul"][x]=(slopesextended[x] * 6) + palt["altjan"][x]
    #altjul has now been fully reconstructed
    
    
    altmonth={}
    for x in range (0,47):
        altmonth["{0}".format(x)]=[] #list order will be 1 to 6 month
    
    #find approx starting altjan points by interpolating for pressure from grid, [jan,jul]
    janinterp=interp1d(palt["pressure"],palt["altjan"])
    julinterp=interp1d(palt["pressure"],palt["altjul"])
    jan2=extrap1d(janinterp)
    jul2=extrap1d(julinterp)
    startingpoints={}
    for x in range (0,47):
        startingpoints["{0}".format(x)]=[]
    for x in range (0,47):
        startingpoints[str(x)].append(int(jan2([pressure[x]])))
        startingpoints[str(x)].append(int(jul2([pressure[x]])))
    
    #now we construct the altitudes for each month based off entries in startingpoints
    for y in range (0,47):
        g_x=interp1d([0,6],startingpoints[str(y)])
        realextrap=extrap1d(g_x)
        for x in range (0,7):
            altmonth[str(y)].append(int(realextrap([x])))
    
    #altmonth now has EVERYTHING we need
    
    files = folders = 0
    for _, dirnames, filenames in os.walk('E:/More netcdf files/'):
        folders += len(dirnames)
    monthlengths=[]
    for x in range (1,folders+1):
        point=str(x)
        if int(x) < 10:
            point = '0' + str(x) #if a number < 10 was entered, put a zero in front of it
        files=0
        for _, dirnames, filenames in os.walk('E:/More netcdf files/' + str(yearselect) + '/'+point+'/'):
            files += len(filenames)
        monthlengths.append(files)
    pyplot.figure(num=yearselect-2003,figsize=(19.2,10.8),dpi=100)
    for y in range (1,folders+1):
        for z in range (1,monthlengths[y-1]+1):
            if y < 10:
                month='0'+str(y)
            else:
                month=str(y)
            if z < 10:
                day='0'+str(z)
            else:
                day=str(z)
            data=Dataset("E:/More netcdf files/" + str(yearselect)+"/" + month + "/ts-"+ str(yearselect)+ month + day + ".nc")
            ocs=data['COS'][0,0:47,0,0]
            for x in range (0,47):
                ocs[x] = ocs[x] * (10 ** 12)
            if y <= 7:
                for x in range (0,47): 
                    pyplot.scatter(ocs[x],altmonth[str(x)][y-1],color=patchcolor)
            if y > 7:
                for x in range (0,47):
                    pyplot.scatter(ocs[x],altmonth[str(x)][5-(y % 8)],color=patchcolor)
    pyplot.legend(handles=[patches.Patch(color=patchcolor,label=str(yearselect))])
    pyplot.xlabel("OCS concentration (ppt)",fontsize=40)
    pyplot.ylabel("Altitude (m)",fontsize=40)
    pyplot.title("OCS Model, " + str(yearselect),fontsize=40)
    axes = pyplot.gca()
    #axes.set_xlim([100,800])
    axes.set_ylim([0,60000])
                    
plotyear(2004)
plotyear(2005)
plotyear(2006)