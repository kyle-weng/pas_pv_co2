#Pasadena and Palos Verdes CO2 Analysis
#Kyle Weng
#July - December 2017

#control flow-----------------------
#sanityCheckAnnualDecomposedScratch
#droughtSplit
#droughtSplitAnnual
#droughtSplitAnnualExtended

### CONTROL GATES
initialWork_Pasadena_PalosVerdes = False
logNormalRun = False
seasonalCycleNonsplitInitiate = False
seasonalCycleNonsplitPlot = False
seasonalCycleSplitInitiate = False
seasonalCycleSplitPlot = False
averagedToOneYear = False
averagedToOneYearDecomposition = False
deviationRatioNonsplitMath = False
deviationRatioNonsplitPlot = False
deviationRatioSplitMath = False
deviationRatioSplitPlot = False
sanityCheckDataHalves = False
sanityCheckLogNormalHalves = False
sanityCheckDecompositionMonthlyData = False
sanityCheckDecompositionMonthlyDataPlot = False
sanityCheckReconstructionMonthlyData = False
sanityCheckReconstructionAspects = False

sanityCheckAnnualDecomposedScratch = True #the entirety of droughtSplit is contingent on this
droughtSplit = True #ie. 2012
droughtSplitPlot = False #dependent on droughtSplit = True
droughtSplitAnnual = True #dependent on droughtSplit
droughtSplitAnnualPlot = False #dependent on droughtSplitAnnual = True
droughtSplitAnnualPlotErrorbars = False #dependent on droughtSplitAnnual = True
droughtSplitAnnualExtended = True #copypaste of dSAP, except extended to 18 months
droughtSplitAnnualTruncated = False #truncated to six months
logNormalDroughtSplit = False
logNormalHistogramsDroughtSplit = False

ml_droughtSetup = False
ml_droughtLogNormalTogether = False
ml_droughtLogNormalHistograms = False

ml_pv_pas_logNormal = False
ml_pv_pas_logNormalSplit = False
lastFifteenYears_dataManipulation = False
lastFifteenYears_histograms = False

CO_O3_correlation = False #this is super time consuming
plot_CO_O3_extraCO2 = False

final_pasadena_three_step = False #dependent on lastFifteenYears_dataManipulation
final_palosverdes_three_step = False

oneYearFit = True #dependent on droughtSplitAnnual

### IMPORTS
from matplotlib import pyplot as plt, patches as pc
from scipy.stats import lognorm,pearsonr
from scipy.integrate import simps
from fractions import Fraction
from calendar import monthrange
from sklearn import linear_model as lm
from datetime import datetime as dt
import statsmodels.api as sm
import numpy as np
import pandas
import time
import math
import sys
import os
import xlsxwriter

global pastime,pvtime

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def ppm(num): #num should be of type numpy.float64 and a scientific notation number (1.6e-05), 
              #turns concentration value to ppm value (*10^-6)
              #returns integer denoting ppm
    xf = Fraction(num)
    return int(round(float(xf.numerator) / float(xf.denominator) * 1000000))

def bestFit(y): #returns the m and b of a line of best fit
                #y = data
                #assume y points are equally spaced from each other in the x direction
    x = np.arange(len(y))
    m,b = np.polyfit(x,y,1)
    return m,b

def nan_helper(y):
    """
    Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def nan_helper_concise(y):
    result = nan_helper(y)
    for x in range(0,len(result)):
        if result == True:
            print x

def interpolateNaN(y):

    nans,x = nan_helper(y)
    y[nans] = np.interp(x(nans),x(~nans),y[~nans])
    return y.round(5)

def removeNaN(*args): #y is a list of datasets from which the same index needs to be removed, y[0] better have the NaN value
    NaNlist = []
    returnList = []
    for x in args:
        for y in range (0,len(x)):
            if pandas.isnull(x[y]) and (not x[y] in NaNlist):
                NaNlist.append(y)
    for x in args:
        for y in range (len(NaNlist)-1,-1,-1):
            x = x[:NaNlist[y]] + x[NaNlist[y]+1:]
        returnList.append(x)
    return returnList
    

def doMath(dset): #interpolate for NaN values, then calculate the best fit and remove the trendline value from the dset
    CO2interpolated = interpolateNaN(dset)
    m,b = bestFit(CO2interpolated)
    CO2removed = removeTrendline(m,b,CO2interpolated)
    return CO2removed

def removeTrendline(m,b,y): #removes value of line of best fit (y = mx + b)
                            #at x point from x index of y data
    z = np.empty(len(y))
    for x in range (0,len(y)):
        z[x] = y[x] - ((m * x) + b)
    return z

def shave(dset,time): #shave first and last NaN values that would be impossible to interpoalte
    
    if (not len(dset) == len(time)):
        print "dset and time need to be the same length"
        return

    counterFront = 0 #number of NaN values (starting from 0) until the first non-NaN element is iterated over
    
    #check front end
    for x in range (0,len(dset)):
        if (not np.isfinite(dset[x])):
            counterFront += 1
        else:
            break
    dset = dset[counterFront:]
    time = time[counterFront:]
    
    #check back end
    counterBack = 0
    for x in range (len(dset)-1,-1,-1):
        if (not np.isfinite(dset[x])):
            counterBack += 1
        else:
            break
    dset = dset[:len(dset)-counterBack]
    time = time[:len(time)-counterBack]
    
    return dset,time

# dset = dataset
# a = start
# b = end
# c = step
# d = name - used for labeling
# e = color
# f = scale of x axis
# g = logspace/linspace
def logNormal(dset,a,b,c,d,e,f,g):
    global shape,loc,scale,distro,pdf
    try:
        plt.xscale(f)
    except ValueError:
        print "Specify the type of scale for the x axis"
        return
    
    if g == "linspace":
        distro = np.linspace(a,b,c)
    elif g == "logspace":
        distro = np.logspace(a,b,c)
        #for x in range (0,len(distro)):
            #print ""
    else:
        print "That didn't work."
        return

    shape,loc,scale = lognorm.fit(dset)
    
    
    pdf = lognorm.pdf(distro, shape, loc, scale)
    
    plt.plot(distro, pdf,color=e)
    
    #INTEGRATE TO CHECK THAT AREA UNDER LOGNORMAL CURVE IS 1
    #Trapezoid
    print "Integration Area Check: " + d
    print "Trapezoidal: " + str(np.trapz(pdf,distro))
    print "Simpson's: " + str(simps(pdf,distro))
    print ""

def logNormal_debug(dset,a,b,c,d,e,f,g): #for debug purposes
    try:
        plt.xscale(f)
    except ValueError:
        print "Specify the type of scale for the x axis"
        return
    
    if g == "linspace":
        distro = np.linspace(a,b,c)
    elif g == "logspace":
        distro = np.logspace(a,b,c)
        #for x in range (0,len(distro)):
            #print ""
    else:
        print "That didn't work."
        return

    shape,loc,scale = lognorm.fit(dset)
    
    
    pdf = lognorm.pdf(distro, shape, loc, scale)
    
    plt.plot(distro, pdf,color=e)
    
    #INTEGRATE TO CHECK THAT AREA UNDER LOGNORMAL CURVE IS 1
    #Trapezoid
    print "Integration Area Check: " + d
    print "Trapezoidal: " + str(np.trapz(pdf,distro))
    print "Simpson's: " + str(simps(pdf,distro))
    print ""
    
    return shape,loc,scale,distro,pdf

def split(dset): #(in half)
    index = int(len(dset) / 2)
    return (dset[:index],dset[index:])

def histogram(dset,bins,color_in):
    #dset should be self explanatory
    #bins is the number of bins you want the data to be sorted into
    #color_in is the color of the histogram
    plt.hist(dset,bins,normed=1,color=color_in) #formerly ax.hist

def pc_regress(x_in,y_in): #IDL code but in Python, not IDL
    global xs,nterm,a,stat
    xs = x_in.shape
    nterm = xs[0]
    a = []
    stat = []
    global tester,results
    
    a = lm.LinearRegression()
    
    a.fit(x_in,y_in)
    
    tester = sm.OLS(y_in,x_in)
    results = tester.fit()
    
    if a.coef_.ndim == 1:
        for iterator in range (0,6):
            if results.bse[iterator] != float(0):
                stat.append(np.abs(a.coef_[iterator]/results.bse[iterator]))
            else:
                stat.append(None)
    elif a.coef_.ndim == 2:
        for iterator in range (0,6):
            if results.bse[iterator] != float(0):
                stat.append(np.abs(a.coef_[0][iterator]/results.bse[iterator]))
            else:
                stat.append(None)
              
    return stat,a,results.bse,a.intercept_,a.coef_
    
def seasonalCycle(dset,time): #IDL code but in Python, not IDL   
    global nt,timeTemp,convertedTimeTemp,convertedTimeDecimal,nn,x,xx,c,r,pcs,y,timeSubtracted
    
    #time conversions
    nt = len(dset)
    
    global convertedTimeTemp #dateTime objects in a list
    
    timeTemp = np.array(time)
    convertedTimeTemp = []
    for x in range (0,len(timeTemp)):
        convertedTimeTemp.append(dt.utcfromtimestamp((timeTemp[x] - np.datetime64('1970-01-01')) / np.timedelta64(1, 's'))) #DeprecationWarning: parsing timezone aware datetimes is deprecated; this will raise an error in the future
    
    global convertedTimeDecimal #floats in a list, ex. 2004.12345
    
    convertedTimeDecimal = []
    for x in range (0,len(convertedTimeTemp)):
        convertedTimeDecimal.append(toYearFraction(convertedTimeTemp[x]))
    
    #difference in time
    nn = (max(convertedTimeDecimal) - min(convertedTimeDecimal))/2
    
    timeSubtracted = [] #ex. subtract 2004 from 2004.1234
    for x in range (0,len(convertedTimeDecimal)):
        timeSubtracted.append(convertedTimeDecimal[x] - convertedTimeDecimal[0])
    
    xx = []
    for x in range (0,len(timeSubtracted)):
        xx.append((timeSubtracted[x]/nn) - 1)
    
    c,r = 6,nt #IDL is column-major
    pcs = np.empty((r,c),dtype='float')
    
    for x in range (0,len(pcs)):
        pcs[x][0] = nn*xx[x]                                #pcs(0,*)
        pcs[x][1] = math.cos(2*math.pi*timeSubtracted[x])   #pcs(1,*)
        pcs[x][2] = math.sin(2*math.pi*timeSubtracted[x])   #pcs(2,*)
        pcs[x][3] = math.cos(4*math.pi*timeSubtracted[x])   #pcs(3,*)
        pcs[x][4] = math.sin(4*math.pi*timeSubtracted[x])   #pcs(4,*)
    
    for x in range (0,nt):
        for y in range (0,len(pcs)):
            pcs[y][5] = 1/3 * nn * nn * 0.5 * (3*(xx[x] ** 2) - 1)
            pcs[y][4] = 1/5 * (nn ** 3) * 0.5 * (5*(xx[x]**3)-3*xx[x])
    
    global coefficients,constant,sigma,trend,season,semi,all_list
    stat, a, sigma, constant, coefficients = pc_regress(pcs,dset)
    if coefficients.ndim == 2:
        coefficients = coefficients[0]
    
    trend = []
    season = []
    semi = []
    all_list = [] #because all is a reserved keyword
    
    for it in range (0,nt):
        trend.append(constant + (coefficients[0] * pcs[it,0]))
        season.append((coefficients[1] * pcs[it,1]) + (coefficients[2] * pcs[it,2]))
        semi.append((coefficients[3] * pcs[it,3]) + (coefficients[4]*pcs[it,4]))
        all_list.append(season[it] + semi[it] + trend[it] + (coefficients[5]*pcs[it,5]))
    
    #The follwoing code corresponds to the above code block
    global amp12,amp6,err12,err6 #these variables aren't used?
    
    amp12 = math.sqrt(coefficients[1] ** 2) + (coefficients[2] ** 2)
    amp6 = math.sqrt(coefficients[3] ** 2) + (coefficients[4] ** 2)
    err12 = math.sqrt(((math.fabs(sigma[1]) + math.fabs(coefficients[1])) ** 2) + ((math.fabs(sigma[2]) + math.fabs(coefficients[2])) ** 2)) - amp12
    err6 = math.sqrt(((math.fabs(sigma[3]) + math.fabs(coefficients[3])) ** 2) + ((math.fabs(sigma[4]) + math.fabs(coefficients[4])) ** 2)) - amp12  
    
    global annamp,semiamp
    
    annamp = []
    semiamp = []
    
    annamp.append(amp12)
    annamp.append(err12)
    
    semiamp.append(amp6)
    semiamp.append(err6)

    return (season,semi,trend,all_list,convertedTimeTemp)

def seasonalCyclePlot(dset,time,title): #instead of calling seasonalCycle() and then plotting, just do it in one go - plots all aspects plus reconstructed
    season,semi,trend,all_l,timeOut = seasonalCycle(dset,time)
    plt.figure()
    plt.plot(timeOut,season,color='blue')
    plt.plot(timeOut,semi,color='red')
    plt.plot(timeOut,trend,color='green')
    plt.plot(timeOut,all_l,color='purple')
    
    #reconstruct
    #all,season,semi,trend
    reconstructed = reconstruct(all_l,semi,season,trend)
    plt.plot(timeOut,reconstructed,color='yellow')
    
    plt.legend(handles=[pc.Patch(label='season',color='blue'),pc.Patch(label='semi',color='red'),pc.Patch(label='trend',color='green'),pc.Patch(label='all',color='purple'),pc.Patch(label='reconstructed',color='yellow')])
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("CO2 (ppm)")
    plt.show()

def seasonalCyclePlot_modified(dset,time,title): #9-13-17 VERSION 2.0 FOR EASY PLOTTING
    def genericAxes2(title):
        plt.title(title)
        plt.xlabel("time")
        plt.ylabel("co2 concentration")
    season,semi,trend,all_l,timeOut = seasonalCycle(dset,time)
    
    plt.figure()
    plt.plot(timeOut,season,color='blue')
    genericAxes2(title + " season")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,semi,color='red')
    genericAxes2(title + " semi")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,trend,color='green')
    genericAxes2(title + " trend")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,all_l,color='purple')
    genericAxes2(title + " all")
    plt.show()
    
    reconstructed = reconstruct(season,semi,trend,all_l)
    decompOriginPlot(timeOut,reconstructed,time,dset,title+" reconstructed")
    
def seasonalCyclePlot_return(dset,time,title): #9-16-17 VERSION 3.0 FOR RETURNING VARIABLES
    def genericAxes2(title):
        plt.title(title)
        plt.xlabel("time")
        plt.ylabel("co2 concentration")
    season,semi,trend,all_l,timeOut = seasonalCycle(dset,time)
    
    plt.figure()
    plt.plot(timeOut,season,color='blue')
    genericAxes2(title + " season")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,semi,color='red')
    genericAxes2(title + " semi")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,trend,color='green')
    genericAxes2(title + " trend")
    plt.show()
    
    plt.figure()
    plt.plot(timeOut,all_l,color='purple')
    genericAxes2(title + " all")
    plt.show()
    
    reconstructed = reconstruct(all_l,season,semi,trend)
    decompOriginPlot(timeOut,reconstructed,time,dset,title+" reconstructed")
    
    return season,semi,trend,all_l,timeOut
    
def genericAxes():
    plt.xlabel('Time (year)')
    plt.ylabel('CO2 concentration (ppm)')
    
def checkConsecutiveDates(df):   #pass dataframe, called by makeConsecutiveDates() so don't call this alone
    missingDates = [] #list of indices that precede missing dates
    for x in range (0,df.shape[0]-1):
        if (df.iloc[x]['date'] + pandas.Timedelta(days=1)) != (df.iloc[x+1]['date']):
            #n is the number of missing days that need to be filled out
            n = (df.iloc[x+1]['date']-df.iloc[x]['date']).days-1
            missingDates.append((x,n))
    if not missingDates: #ie. if it is empty
        return (True,missingDates)
    else:
        return (False,missingDates)
        
def makeConsecutiveDates(df_in):   #pass dataframe
    #global toBeInserted,df_modified,result,index,repeat,toBeInsertedList,x,z
    
    df_modified = pandas.DataFrame()
    
    result = checkConsecutiveDates(df_in)
    if result[0] == True:
        print "Dates are already consecutive!"
        return df_in
    toBeInsertedList = []
    for x in range (0,len(result[1])):
        index = result[1][x][0]
        repeat = result[1][x][1]
        toBeInserted = pandas.DataFrame()
        for y in range (1,repeat+1):
             toBeInserted = pandas.concat([toBeInserted,pandas.DataFrame({'date':[df_in.iloc[index]['date'] + pandas.Timedelta(days=y)],'co2':None})])
        toBeInsertedList.append(toBeInserted)
    
    #append method
    z = 0 #index of toBeInsertedList
    y = 0 #index of original dataframe, df_in
    for x in range (0,len(toBeInsertedList)+df_in.shape[0]):
        ab = df_in.iloc[y:y+1]
        df_modified = df_modified.append(ab)
        y += 1
        if z < len(toBeInsertedList):
            if x == result[1][z][0]:
                insert = toBeInsertedList[z]
                z += 1
                df_modified = df_modified.append(insert)
    df_modified.reset_index(drop=True)
    
    df_reindexed = pandas.DataFrame({'date':np.array(df_modified['date']),'co2':np.array(df_modified['co2'])},columns=['date','co2'])
    
    return df_reindexed

def dayOfYear(timestamp):
    x = dt.strptime(str(timestamp),"%Y-%m-%d %H:%M:%S")
    return x.timetuple().tm_yday

def averageYears_day(df_in): #returns one year of averaged year data
    global df,month_day
    co2_list = []
    month_day = []
    
    #THROWING AWAY DATA FROM 2/29 ON LEAP YEARS
    
    for x in range (1,366):
        co2_list.append([])
    
    for x in range (1,13):
        for y in range (1,monthrange(2001,x)[1]+1):
            string_to_append = ""
            if x < 10:
                string_to_append += "0"
            string_to_append += (str(x) + "_")
            if y < 10:
                string_to_append += "0"
            string_to_append += str(y)
            month_day.append(string_to_append)
    
    df = pandas.DataFrame({'date':month_day,'co2 list':co2_list},columns=['date','co2 list'])
    
    for x in range (0,len(df_in)):
        if df_in.iloc[x]['date'].month == 2 and df_in.iloc[x]['date'].day == 29:
            continue
        else:
            if df_in.iloc[x]['date'].year % 4 == 0 and (df_in.iloc[x]['date'] > dt(df_in.iloc[x]['date'].year,2,29)):
                indx = dayOfYear(df_in.iloc[x]['date']) - 2
            else:
                indx = dayOfYear(df_in.iloc[x]['date']) - 1
        df.iloc[indx]['co2 list'].append(df_in.iloc[x]['co2'])
    
    averaged_by_day = []
    for x in range (0,len(df)):
        averaged_by_day.append(np.mean(df.iloc[x]['co2 list']))
    
    return averaged_by_day

def averageYears_month(df_in):
    monthData = []
    monthlyAverages = []
    for x in range (1,13):
        monthData.append([])
    for x in range (0,len(df_in)):
        monthData[df_in.iloc[x]['date'].month-1].append(df_in.iloc[x]['co2'])
    for x in range (1,13):
        monthlyAverages.append(np.mean(monthData[x-1]))
    return monthlyAverages

def seasonalLegend():
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena'),pc.Patch(color = 'red',label = 'Palos Verdes')])

def decomposePlot(pasTimeProper,pasDecompose,pvTimeProper,pvDecompose,title): #plot decomposed data, Pasadena and PV together
                #Ex. decomposePlot(<pas time>,<pas decomposed data>,<pv time>,<pv decomposed       
                #data>,<plot title>)
    plt.figure()
    plt.plot(pasTimeProper,pasDecompose,'blue')
    plt.plot(pvTimeProper,pvDecompose,'red')
    seasonalLegend()
    plt.title(title)
    genericAxes()
    plt.show()

def removeBefore2004(dset,time):
    for x in range (0,len(time)):
        if str(time[x])[:4] == '2004':
            y = x
            break
    return (dset[y:],time[y:])

def splitBeforeYear(dset,time,year):
    for x in range (0,len(time)):
        if str(time.iloc[x])[:4] == str(year):
            y = x
            break
    return ((dset[y:],time[y:]),(dset[:y],time[:y]))

def deviation_ratio(co2DSet): #dCO2,rCO2
    interpd = interpolateNaN(co2DSet)
    interpd = np.array(interpd)
    dSetMinusAverage = []
    dCO2 = []
    rCO2 = []
    dSetAverage = np.mean(interpd)
    for x in range (0,len(interpd)):
        dSetMinusAverage.append(interpd[x]-dSetAverage)
    minDset = min(dSetMinusAverage)
    for x in range (0,len(interpd)):
        dCO2.append(interpd[x]-dSetAverage-minDset+1e-20)
        rCO2.append(interpd[x]/dSetAverage)
    return interpd,dSetMinusAverage,dSetAverage,dCO2,rCO2

def monthlyAverageTimeSetup():
    
    PVMonthlyTime = []
    PasMonthlyTime = []
    
    for x in range (6,14):
        if x == 6:
            limLower = 2
        else:
            limLower = 1
        for y in range (limLower,13):
            toBeAddedMonth = str(y)
            toBeAddedYear = str(x)
            
            if y < 10:
                toBeAddedMonth = "0" + toBeAddedMonth
            if x < 10:
                toBeAddedYear = "0" + toBeAddedYear
            
            PasMonthlyTime.append(pandas._libs.tslib.Timestamp(toBeAddedMonth + "-01-20" + toBeAddedYear))
            
    for x in range (9,14):
        if x == 9:
            limLower = 10
        else:
            limLower = 1
        for y in range (limLower,13):
            toBeAddedMonth = str(y)
            toBeAddedYear = str(x)
            
            if y < 10:
                toBeAddedMonth = "0" + toBeAddedMonth
            if x < 10:
                toBeAddedYear = "0" + toBeAddedYear
            
            PVMonthlyTime.append(pandas._libs.tslib.Timestamp(toBeAddedMonth + "-01-20" + toBeAddedYear))
            
    return (PasMonthlyTime,PVMonthlyTime)

def genericPlot(x,y,title):
    plt.figure()
    plt.plot(x,y)
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("CO2 ppm")
    plt.show()

def reconstruct(dAll,dSeason,dSemi,dTrend):
    reconstructed = []
    if not(len(dAll) == len(dSeason) == len(dSemi) == len(dTrend)):
        print "reconstruction didn't work; all inputs need to be the same length"
        return
    for x in range (0,len(dAll)):
        reconstructed.append(dSeason[x] + dSemi[x] + dTrend[x])
    return reconstructed

def decompOriginPlot(decomposedTime,decomposedReconstructed,actualTime,actualData,title):
    plt.figure()
    plt.plot(actualTime,actualData,color = 'orange')
    plt.plot(decomposedTime,decomposedReconstructed,'--',color = 'blue')
    plt.title(title)
    plt.legend(handles=[pc.Patch(label='decomposed',linestyle='--',color='blue'),pc.Patch(label='actual',color = 'orange')])
    genericAxes()
    plt.xticks(list(np.arange(1,len(decomposedTime)+1,1)))
    plt.show()

def decompOriginPlotErrorbars(decomposedTime,decomposedReconstructed,actualTime,actualData,title,error):
    plt.figure()
    plt.errorbar(actualTime,actualData,yerr=error,color = 'orange',fmt='o',capsize=3)
    plt.plot(decomposedTime,decomposedReconstructed,'--',color = 'blue')
    plt.title(title)
    plt.legend(handles=[pc.Patch(label='decomposed',linestyle='--',color='blue'),pc.Patch(label='actual',color = 'orange')])
    genericAxes()
    plt.xticks(list(np.arange(1,len(decomposedTime)+1,1)))
    plt.show()

def makeDTYear(year,TS_type):
    dty = []
    for x in range (1,13):
        for y in range (1,monthrange(year,x)[1]+1):
            dty.append(dt(year,x,y))
    if TS_type == "numpy":
        return dty
    elif TS_type == "pandas":
        return pandas.core.series.Series(dty)
    else:
        print "need to specify numpy or pandas timestamps"
        return

def makeDTYear_m(year,TS_type):
    dty=[]
    for x in range (1,13):
        dty.append(dt(year,x,1))
    if TS_type == "numpy":
        return dty
    elif TS_type == "pandas":
        return pandas.core.series.Series(dty)
    else:
        print "need to specify numpy or pandas timestamps"
        return

def plot1dSeries(inpt): #debugging purposes
    plt.figure()
    plt.plot(np.array(range(0,len(inpt))),inpt)
    plt.show()

def averageYearMonth(dset,time,year,month):
    global byYearMonthIndices,byYearMonth
    if str(type(time)) == "<class 'pandas.core.series.Series'>":
        byYearMonthIndices = []
        byYearMonth = []
        for x in range (0,len(time)):
            if time.iloc[x].month == month and time.iloc[x].year == year:
                byYearMonthIndices.append(x)
        for x in range (0,len(byYearMonthIndices)):
            byYearMonth.append(dset[byYearMonthIndices[x]])
        return byYearMonth
    if str(type(time)) == "<type 'list'>":
        byYearMonthIndices = []
        byYearMonth = []
        for x in range (0,len(time)):
            if time[x].month == month and time[x].year == year:
                byYearMonthIndices.append(x)
        for x in range (0,len(byYearMonthIndices)):
            byYearMonth.append(dset[byYearMonthIndices[x]])
        return byYearMonth
def averageMonth(dset,time,month):
    byMonthIndices = []
    byMonth = []
    for x in range (0,len(time)):
        if time.iloc[x].month == month:
            byMonthIndices.append(x)
    for x in range (0,len(byMonthIndices)):
        byMonth.append(dset[x])
    return byMonth

def annualPlotMonthPoint(dset,timeData,month): #9-18-17
    global monthDataAllYears
    if str(type(timeData)) == "<class 'pandas.core.series.Series'>":
        #dset = dset.reset_index(drop=True)
        timeData = timeData.reset_index(drop=True)
        monthDataAllYears = [] #0 index is actual data, 1 index is mean of said data
        lowerYearBound = timeData.iloc[0].year
        upperYearBound = timeData.iloc[timeData.size-1].year
        for x in range (lowerYearBound,upperYearBound+1):
            monthDataAllYears.append((averageYearMonth(dset,timeData,x,month),np.mean(averageYearMonth(dset,timeData,x,month))))
        #shave off any end bits for empty month sets (ex. in the march dset if there is no march 2017 data)
        for x in range (len(monthDataAllYears)-1,-1,-1):
            if np.isnan(monthDataAllYears[x][1]):
                monthDataAllYears.remove(monthDataAllYears[x])
        return np.mean([x[1] for x in monthDataAllYears]),[np.mean(x[0]) for x in monthDataAllYears]
    if str(type(timeData)) == "<type 'list'>":
        monthDataAllYears = [] #0 index is actual data, 1 index is mean of said data
        lowerYearBound = timeData[0].year
        upperYearBound = timeData[len(timeData)-1].year
        for x in range (lowerYearBound,upperYearBound+1):
            monthDataAllYears.append((averageYearMonth(dset,timeData,x,month),np.mean(averageYearMonth(dset,timeData,x,month))))
        for x in range (len(monthDataAllYears)-1,-1,-1):
            if np.isnan(monthDataAllYears[x][1]):
                monthDataAllYears.remove(monthDataAllYears[x])
        #return np.mean([x[1] for x in monthDataAllYears]),[np.mean(x[0]) for x in monthDataAllYears]
        return [np.mean(x[0]) for x in monthDataAllYears]

def repeater(dset,desired_length): #truncates or repeats LISTS as necessary
    if desired_length > len(dset):
        result_to_return = []
        for x in range (0,desired_length / len(dset)):
            result_to_return += dset
        result_to_return += dset[:desired_length%len(dset)]
        return result_to_return
    elif desired_length == len(dset):
        return dset
    elif desired_length < len(dset):
        return dset[:desired_length]

def removeDataPoint(dset,dset_time,point_time_to_remove):
    #dset = Pandas object
    #dset_time = Pandas object, must match dset in length
    #point_time_to_remove must be of type pandas.Timestamp() (ex. pandas.Timestamp('2014-08-30 00:00:00')) in dset_time
    dset_time = dset_time.reset_index(drop = True)
    
    if len(dset) != len(dset_time):
        print "removeDataPoint() didn't work because of unequal input array length"
        return
    
    if str(type(dset)) == "<type 'numpy.ndarray'>":
        return_dset = np.delete(dset,dset_time.loc[dset_time==point_time_to_remove].index.values[0])
    elif str(type(dset)) == "<class 'pandas.core.series.Series'>":
        dset = dset.reset_index(drop = True)
        return_dset = dset.drop(dset_time.loc[dset_time==point_time_to_remove].index.values[0])
    
    if str(type(dset_time)) == "<class 'pandas.core.series.Series'>":
        return_dset_time = dset_time.drop(dset_time.loc[dset_time==point_time_to_remove].index.values[0])
    else:
        print "dset_time needs to be a pandas series"
        return
    
    return (return_dset,return_dset_time)

def nanValueShave_end(dset): #SHAVE O3 AND CO DATASETS AT THE END OF THE ACTUAL DATA
    for x in range (0,len(dset)):
        if np.isnan(dset.iloc[x]["value"]):
            break
    return dset[:x]

def findAdjacentDates(date_to_place,date_data): #find dates adjacent to date_to_place in date_data
    #date_to_place: single datetime or timestamp object
    #date_data: datetime or timestamp list/array, this better be sorted low/high
    for x in range (0,len(date_data)-1):
        if date_data[x+1] > date_to_place and date_data[x] < date_to_place:
            return date_data[x:x+2]
        elif date_data[x] == date_to_place:
            return date_data[x]
        if x == len(date_data)-1:
            if date_data[x+1] == date_to_place:
                return date_data[x+1]
    return "findAdjacentDates failed"

def listToExcel(filename,*args):
    if type("hello") != str:
        print "listToExcel needs a string filename"
        return
    
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()
    worksheet.set_column()
    #THE ABOVE IS NOT FINISHED

#def nanValueSubtractionShave(dset): #PANDAS DATAFRAMES ONLY, CO AND O3 ONLY
    

data = pandas.ExcelFile("C:\Users\Kyle W\Desktop\Pas PV flasks CO2 13C 18O summ for Kyle.xlsx") #local
#data = pandas.ExcelFile("/home/kaw/Pasadena_PalosVerdes_CO2/Pas PV flasks CO2 13C 18O summ for Kyle.xlsx") #server

pasadena_flasks = data.parse("Pasadena_flasks")
pasCO2 = pasadena_flasks["CO2 (ppm)"]
pasTime = pasadena_flasks["date"]

palosverdes_flasks = data.parse("PV_flasks_aves")
pvCO2 = palosverdes_flasks["[CO2]B"]
pvTime = palosverdes_flasks["Date"]

"""
laJollaData = pandas.ExcelFile("C:\Users\Kyle W\Desktop\dailyflasks_co2_ljo 07302017.xlsx")
lajolla_flasks = laJollaData.parse("data")
ljCO2 = lajolla_flasks["good CO2 (ppm)"]
ljDate = lajolla_flasks["Excel date"]
"""

laJollaData = pandas.read_csv("C:\Users\Kyle W\Desktop\daily_flask_co2_ljo.csv",skiprows=67)

maunaLoaData = pandas.ExcelFile("C:\Users\Kyle W\Desktop\MLO CO2 1969-2016 flasks.xlsx")
maunaloa_flasks = maunaLoaData.parse("data")
#BRIEF PROCESSING - GET RID OF ALL ROWS WITH CO2 NAN MEASUREMENT
maunaloa_flasks = maunaloa_flasks[np.isfinite(maunaloa_flasks['good CO2 (ppm)'])]
mlCO2 = maunaloa_flasks["good CO2 (ppm)"]
mlDate = maunaloa_flasks["date and time"]

### Pasadena and Palos Verdes-- initial math: Shave NaNs, interpolate, calculate + remove trendline
if initialWork_Pasadena_PalosVerdes:
    pasCO2Shaved,pasTimeShaved = shave(pasCO2,pasTime)
    pvCO2Shaved,pvTimeShaved = shave(pvCO2,pvTime)
    
    #interpolate for NaN values, then calculate the best fit and remove the trendline value from the dset
    
    pasCO2removed = doMath(pasCO2Shaved)
    pvCO2removed = doMath(pvCO2Shaved)
    """
    pasCO2removed = interpolateNaN(pasCO2Shaved)
    pvCO2removed = interpolateNaN(pvCO2Shaved)
    
    for x in range (0,len(pasCO2removed)):
        pasCO2removed[x] -= np.mean(pasCO2removed)
    for x in range (0,len(pvCO2removed)):
        pvCO2removed[x] -= np.mean(pvCO2removed)
    """
    
    pas2004CO2,pas2004Time = removeBefore2004(pasCO2removed,pasTimeShaved)
    pasFirstTime,pasSecondTime = split(pas2004Time)
    pasFirstHalf,pasSecondHalf = split(pas2004CO2)
    pvFirstHalf,pvSecondHalf = split(pvCO2removed)

### logNormal
if logNormalRun:
    plt.figure()
    logNormal(pasFirstHalf,0,3,1000,"Pasadena First Half",'blue','log','logspace')
    logNormal(pasSecondHalf,0,3,1000,"Pasadena Second Half",'red','log','logspace')
    logNormal(pvFirstHalf,0,3,1000,"PV First Half",'green','log','logspace')
    logNormal(pvSecondHalf,0,3,1000,"PV Second Half",'orange','log','logspace')
    plt.title("PV + Pas log normal, log scale")
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena 1st Half'),pc.Patch(color='red',label='Pasadena 2nd Half'),pc.Patch(color = 'green',label = 'PV 1st Half'),pc.Patch(color = 'orange',label = 'PV 2nd Half')])
    plt.show()
    
    plt.figure()
    logNormal(pasFirstHalf,-40,120,100,"Pasadena First Half",'blue','linear','linspace')
    logNormal(pasSecondHalf,-40,120,100,"Pasadena Second Half",'red','linear','linspace')
    logNormal(pvFirstHalf,-40,120,100,"PV First Half",'green','linear','linspace')
    logNormal(pvSecondHalf,-40,120,100,"PV Second Half",'orange','linear','linspace')
    plt.title("PV + Pas log normal, lin scale")
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena 1st Half'),pc.Patch(color='red',label='Pasadena 2nd Half'),pc.Patch(color = 'green',label = 'PV 1st Half'),pc.Patch(color = 'orange',label = 'PV 2nd Half')])
    plt.show()

### NONSPLIT SEASONAL CYCLES
if seasonalCycleNonsplitInitiate:
    pasSeason,pasSemi,pasTrend,pasAll,pasTimeProper = seasonalCycle(pas2004CO2,pas2004Time)
    pvSeason,pvSemi,pvTrend,pvAll,pvTimeProper = seasonalCycle(pvCO2removed,pvTime)
if seasonalCycleNonsplitPlot:
    #season
    decomposePlot(pasTimeProper,pasSeason,pvTimeProper,pvSeason,"PV/Pas Season")
    
    #semi
    decomposePlot(pasTimeProper,pasSemi,pvTimeProper,pvSemi,"PV/Pas Semi")
    
    #trend
    decomposePlot(pasTimeProper,pasTrend,pvTimeProper,pvTrend,"PV/Pas Trend")
    
    #all-- this is special so i'm not sticking it in decomposePlot()
    plt.figure()
    plt.plot(pasTimeProper,pasAll,'blue')
    plt.plot(pvTimeProper,pvAll,'red')
    #seasonalLegend()
    plt.title("PV/Pas All")
    genericAxes()
    plt.show()
    #add the original data as well
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena Decomposed'),pc.Patch(color = 'red',label = 'Palos Verdes Decomposed'),pc.Patch(color = 'green',label = 'Pasadena Data'),pc.Patch(color='purple',label='Palos Verdes Data')])
    
    
    #All, but separate for PV and Pas (plot all and original data)
    plt.figure()
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pasTimeProper,pasAll,'blue')
    plt.legend(handles=[pc.Patch(color='blue',label='Pasadena Decomposed'),pc.Patch(color='green',label='Pasadena Data')])
    genericAxes()
    plt.title("Pasadena, Decomposed vs. Data")
    plt.show()
    
    plt.figure()
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.plot(pvTimeProper,pvAll,'red')
    plt.legend(handles=[pc.Patch(color='red',label='Palos Verdes Decomposed'),pc.Patch(color='purple',label='Palos Verdes Data')])
    genericAxes()
    plt.title("Palos Verdes, Decomposed vs. Data")
    plt.show()


### SPLIT SEASONAL CYCLES
if seasonalCycleSplitInitiate:
    #Split first
    pasSplitTime = split(pas2004Time)
    pasSplitCO2 = split(pas2004CO2)
    
    pvSplitTime = split(pvTime)
    pvSplitCO2 = split(pvCO2removed)
    
    #remove mean 1.5th
    
    #Seasonal cycle second
    pasSeason1,pasSemi1,pasTrend1,pasAll1,pasTimeProper1 = seasonalCycle(pasSplitCO2[0],pasSplitTime[0])
    pasSeason2,pasSemi2,pasTrend2,pasAll2,pasTimeProper2 = seasonalCycle(pasSplitCO2[1],pasSplitTime[1])
    
    pvSeason1,pvSemi1,pvTrend1,pvAll1,pvTimeProper1 = seasonalCycle(pvSplitCO2[0],pvSplitTime[0])
    pvSeason2,pvSemi2,pvTrend2,pvAll2,pvTimeProper2 = seasonalCycle(pvSplitCO2[1],pvSplitTime[1])

if seasonalCycleSplitPlot:
    #season
    decomposePlot(pasTimeProper1,pasSeason1,pvTimeProper1,pvSeason1,"PV/Pas Season First Half")
    decomposePlot(pasTimeProper2,pasSeason2,pvTimeProper2,pvSeason2,"PV/Pas Season Second Half")
    
    #semi
    decomposePlot(pasTimeProper1,pasSemi1,pvTimeProper1,pvSemi1,"PV/Pas Semi First Half")
    decomposePlot(pasTimeProper2,pasSemi2,pvTimeProper2,pvSemi2,"PV/Pas Semi Second Half")
    
    #trend
    decomposePlot(pasTimeProper1,pasTrend1,pvTimeProper1,pvTrend1,"PV/Pas Trend First Half")
    decomposePlot(pasTimeProper2,pasTrend2,pvTimeProper2,pvTrend2,"PV/Pas Trend Second Half")
    
    
                # All-- data batch 1
                    # all
    
    plt.figure()
    plt.plot(pasTimeProper1,pasAll1,'blue')
    plt.plot(pvTimeProper1,pvAll1,'red')
    seasonalLegend()
    plt.title("PV/Pas All First Half")
    genericAxes()
    plt.show()
    """
                        ### add the original data as well
    
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena Decomposed'),pc.Patch(color = 'red',label = 'Palos Verdes Decomposed'),pc.Patch(color = 'green',label = 'Pasadena Data'),pc.Patch(color='purple',label='Palos Verdes Data')])
    
                # All, but separate for PV and Pas (plot all and original data)
    
    plt.figure()
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pasTimeProper1,pasAll1,'blue')
    plt.legend(handles=[pc.Patch(color='blue',label='Pasadena Decomposed'),pc.Patch(color='green',label='Pasadena Data')])
    genericAxes()
    plt.title("Pasadena, Decomposed vs. Data")
    plt.show()
    
    plt.figure()
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.plot(pvTimeProper,pvAll,'red')
    plt.legend(handles=[pc.Patch(color='red',label='Palos Verdes Decomposed'),pc.Patch(color='purple',label='Palos Verdes Data')])
    genericAxes()
    plt.title("Palos Verdes, Decomposed vs. Data")
    plt.show()
    """
    
    
                # All-- data batch 2
                    # all
    
    plt.figure()
    plt.plot(pasTimeProper2,pasAll2,'blue')
    plt.plot(pvTimeProper2,pvAll2,'red')
    seasonalLegend()
    plt.title("PV/Pas All Second Half")
    genericAxes()
    plt.show()
    """
    #add the original data as well
    
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena Decomposed'),pc.Patch(color = 'red',label = 'Palos Verdes Decomposed'),pc.Patch(color = 'green',label = 'Pasadena Data'),pc.Patch(color='purple',label='Palos Verdes Data')])
    
                # All, but separate for PV and Pas (plot all and original data)
    
    plt.figure()
    plt.plot(pas2004Time,pas2004CO2,'green')
    plt.plot(pasTimeProper,pasAll,'blue')
    plt.legend(handles=[pc.Patch(color='blue',label='Pasadena Decomposed'),pc.Patch(color='green',label='Pasadena Data')])
    genericAxes()
    plt.title("Pasadena, Decomposed vs. Data")
    plt.show()
    
    plt.figure()
    plt.plot(pvTime,pvCO2removed,'purple')
    plt.plot(pvTimeProper,pvAll,'red')
    plt.legend(handles=[pc.Patch(color='red',label='Palos Verdes Decomposed'),pc.Patch(color='purple',label='Palos Verdes Data')])
    genericAxes()
    plt.title("Palos Verdes, Decomposed vs. Data")
    plt.show()
    """


### PV and Pas, averaged to one year (requires seasonalCycleNonsplitInitiate)
if averagedToOneYear:
    dFPasDecomposed = pandas.DataFrame({'date':pasTimeProper,'co2':pasAll},columns=['date','co2'])
    dFPasActual = pandas.DataFrame({'date':list(pas2004Time),'co2':list(pas2004CO2)},columns=['date','co2'])
    dFrame2 = makeConsecutiveDates(dFPasDecomposed)
    dFrame3 = makeConsecutiveDates(dFPasActual)    
    
    dFrame2['co2'] = interpolateNaN(dFrame2['co2'])
    dFrame3['co2'] = interpolateNaN(dFrame3['co2'])
    
    pasDataAveraged = averageYears_day(dFrame3)
    pasAllAveraged = averageYears_day(dFrame2) #decomposed
    
    dFPVDecomposed = pandas.DataFrame({'date':pvTimeProper,'co2':pvAll},columns=['date','co2'])
    dFPVActual =pandas.DataFrame({'date':list(pvTime),'co2':list(pvCO2removed)},columns=['date','co2'])
    
    dFrame4 = makeConsecutiveDates(dFPVDecomposed)
    dFrame5 = makeConsecutiveDates(dFPVActual)
    dFrame4['co2'] = interpolateNaN(dFrame4['co2'])
    dFrame5['co2'] = interpolateNaN(dFrame5['co2'])
    
    pvDataAveraged = averageYears_day(dFrame5)
    pvAllAveraged = averageYears_day(dFrame4) #decomposed
    
    handles_in=[pc.Patch(color='blue',label='Data'),pc.Patch(color='green',label='Decomposed')]
    
    plt.figure()
    plt.title("Pas CO2, averaged to one year")
    plt.plot(np.arange(1,366),pasDataAveraged,color='blue')
    plt.plot(np.arange(1,366),pasAllAveraged,color='green')
    plt.legend(handles=handles_in)
    plt.ylabel("CO2 (ppm)")
    plt.xlabel("Day of year (Jan 1 = 0)")
    plt.show()
    
    plt.figure()
    plt.title("PV CO2, averaged to one year")
    plt.plot(np.arange(1,366),pvDataAveraged,color='blue')
    plt.plot(np.arange(1,366),pvAllAveraged,color='green')
    plt.legend(handles=handles_in)
    plt.ylabel("CO2 (ppm)")
    plt.xlabel("Day of year (Jan 1 = 0)")
    plt.show()

### PV and Pas, averaged to one year decomposition
if averagedToOneYearDecomposition:
    yearSelector = 2001 #literally just put any year here
    yearDT = makeDTYear(yearSelector,"numpy")
    #pvSeason1,pvSemi1,pvTrend1,pvAll1,pvTimeProper1 = seasonalCycle(pvSplitCO2[0],pvSplitTime[0])
    pasYearSeason,pasYearSemi,pasYearTrend,pasYearAll,pasYearTime = seasonalCycle(pasDataAveraged,pandas.Series(yearDT))
    pvYearSeason,pvYearSemi,pvYearTrend,pvYearAll,pvYearTime = seasonalCycle(pvDataAveraged,pandas.Series(yearDT))
    genericPlot(pasYearTime,pasYearSeason,"pas yearly season")
    genericPlot(pasYearTime,pasYearSemi,"pas yearly semi")
    genericPlot(pasYearTime,pasYearTrend,"pas yearly trend")
    genericPlot(pasYearTime,pasYearAll,"pas yearly all")
    
    genericPlot(pvYearTime,pvYearSeason,"pv yearly season")
    genericPlot(pvYearTime,pvYearSemi,"pv yearly semi")
    genericPlot(pvYearTime,pvYearTrend,"pv yearly trend")
    genericPlot(pvYearTime,pvYearAll,"pv yearly all")

#dCO2, rCO2
""" From these CO2 mixing ratio (ppm) data, calculate log(CO2), log(dCO2) and log(rCO2) and their distributions. Find out which one is closest to the normal distribution. dCO2 = CO2 - average(CO2) - min[CO2 - average(CO2)] + 1.e-20 to make sure dCO2>0 and rCO2 = CO2/average(CO2) """

### dCO2, rCO2 NONSPLIT
if deviationRatioNonsplitMath:
    pasCO2interpolated,pasDataMinusAverage,pasDataAverage,dCO2Pas,rCO2Pas = deviation_ratio(pasCO2)
    pvCO2interpolated,pvDataMinusAverage,pvDatatAverage,dCO2PV,rCO2PV = deviation_ratio(pvCO2)

### dCO2, rCO2 SPLIT
if deviationRatioSplitMath:
    pvCO2FirstHalf,pvCO2SecondHalf = split(pvCO2)
    pasCO2FirstHalf,pasCO2SecondHalf = split(pasCO2)
    
    #First Half = FH
    pvCO2interpolated_FH,pvDataMinusAverage_FH,pvDatatAverage_FH,dCO2PV_FH,rCO2PV_FH = deviation_ratio(pvCO2FirstHalf)
    pasCO2interpolated_FH,pasDataMinusAverage_FH,pasDataAverage_FH,dCO2Pas_FH,rCO2Pas_FH = deviation_ratio(pasCO2FirstHalf)
    
    #Second Half = SH
    pvCO2interpolated_SH,pvDataMinusAverage_SH,pvDatatAverage_SH,dCO2PV_SH,rCO2PV_SH = deviation_ratio(pvCO2SecondHalf)
    pasCO2interpolated_SH,pasDataMinusAverage_SH,pasDataAverage_SH,dCO2Pas_SH,rCO2Pas_SH = deviation_ratio(pasCO2SecondHalf)

### dCO2 and rCO2 by d/r (between sites), nonsplit data
if deviationRatioNonsplitPlot:
        #dCO2
    plt.figure()
    logNormal(dCO2PV,-1,2.5,1000,"dCO2",'red','log','logspace')
    logNormal(dCO2Pas,-1,2.5,1000,"dCO2",'blue','log','logspace')
    plt.legend(handles=[pc.Patch(color='red',label='Palos Verdes'),pc.Patch(color='blue',label='Pasadena')])
    plt.title("dCO2 between sites")
        #rCO2
    plt.figure()
    logNormal(rCO2PV,0.75,1.25,1000,"dCO2",'red','linear','linspace')
    logNormal(rCO2Pas,0.75,1.25,1000,"dCO2",'blue','linear','linspace')
    #plt.axvline(x=1)
    plt.legend(handles=[pc.Patch(color='red',label='Palos Verdes'),pc.Patch(color='blue',label='Pasadena')])
    plt.title("rCO2 between sites")

### dCO2 and rCO2 by d/r (between sites), split data
if deviationRatioSplitPlot:
        #dCO2, zoomed out for dCO2PV_FH
    handles_in = []
    plt.figure()
    
    logNormal(dCO2PV_FH,-4,2.5,1000,"dCO2PV_FH",'red','log','logspace')
    handles_in.append(pc.Patch(color='red',label='Palos Verdes First Half'))
    
    logNormal(dCO2PV_SH,-4,2.5,1000,"dCO2PV_SH",'blue','log','logspace')
    handles_in.append(pc.Patch(color='blue',label = 'Palos Verdes Second Half'))
    
    logNormal(dCO2Pas_FH,-4,2.5,1000,"dCO2Pas_FH",'green','log','logspace')
    handles_in.append(pc.Patch(color='green',label='Pasadena First Half'))
    
    logNormal(dCO2Pas_SH,-4,2.5,1000,"dCO2Pas_SH",'purple','log','logspace')
    handles_in.append(pc.Patch(color='purple',label='Pasadena Second Half'))
    plt.legend(handles=handles_in)
    plt.title("dCO2 between sites, split by halves, zoomed out")
    
        #dCO2, zoomed in (normal view)
    plt.figure()
    
    logNormal(dCO2PV_FH,-0.05,2.3,1000,"dCO2PV_FH",'red','log','logspace')
    logNormal(dCO2PV_SH,-0.05,2.3,1000,"dCO2PV_SH",'blue','log','logspace')
    logNormal(dCO2Pas_FH,-0.05,2.3,1000,"dCO2Pas_FH",'green','log','logspace')
    logNormal(dCO2Pas_SH,-0.05,2.3,1000,"dCO2Pas_SH",'purple','log','logspace')
    plt.legend(handles=handles_in)
    plt.title("dCO2 between sites, split by halves, zoomed in")
    
        #rCO2
    plt.figure()
    
    logNormal(rCO2PV_FH,0.87,1.2,1000,"rCO2PV_FH",'red','linear','linspace')
    logNormal(rCO2PV_SH,0.87,1.2,1000,"rCO2PV_SH",'blue','linear','linspace')
    logNormal(rCO2Pas_FH,0.87,1.2,1000,"rCO2Pas_FH",'green','linear','linspace')
    logNormal(rCO2Pas_SH,0.87,1.2,1000,"rCO2Pas_SH",'purple','linear','linspace')
    plt.legend(handles=handles_in)
    plt.title("rCO2 between sites, split by halves")

### sanity check-- halves of data
if sanityCheckDataHalves:
    plt.figure()
    plt.plot(pvCO2interpolated_FH,color='blue')
    plt.plot(pvCO2interpolated_SH,color='red')
    plt.title("PV")
    plt.legend(handles=[pc.Patch(color='blue',label='first half'),pc.Patch(color='red',label='second half')])
    plt.show()
    
    plt.figure()
    plt.plot(pasCO2interpolated_FH,color='blue')
    plt.plot(pasCO2interpolated_SH,color='red')
    plt.legend(handles=[pc.Patch(color='blue',label='first half'),pc.Patch(color='red',label='second half')])
    plt.title("Pas")
    plt.show()

### sanity check-- logNormal_debug halves
if sanityCheckLogNormalHalves:
    PVFH_shape,PVFH_loc,PVFH_scale,PVFH_distro,PVFH_pdf = logNormal_debug(dCO2PV_FH,-4,2.5,1000,"dCO2PV_FH",'red','log','logspace')
    PVSH_shape,PVSH_loc,PVSH_scale,PVSH_distro,PVSH_pdf = logNormal_debug(dCO2PV_SH,-4,2.5,1000,"dCO2PV_SH",'red','log','logspace')
    
### sanity check-- decomposition of monthly data
if sanityCheckDecompositionMonthlyData:
    PVMonthly = pandas.read_table("C:\Users\Kyle W\Desktop\CO2_PAS_PV\PV_co2ff_2009_13_06252015.txt")
    PasMonthly = pandas.read_table("C:\Users\Kyle W\Desktop\CO2_PAS_PV\Pas_co2ff06252015.txt")
    
    PasMonthlyTime,PVMonthlyTime = monthlyAverageTimeSetup()
    
    PVMonthlyAvg = np.mean(PVMonthly)
    for x in range (0,len(PVMonthly)):
        PVMonthly.iloc[x]["CO2ff: Monthly data are from Oct 2009 through Dec 2013"] -= PVMonthlyAvg
    
    PasMonthlyAvg = np.mean(PasMonthly)
    for x in range (0,len(PasMonthly)):
        PasMonthly.iloc[x]["[CO2]ff deseason: Monthly average data are from Feb 19, 2006 to Dec 2013"] -= PasMonthlyAvg
    
    pasMSeason,pasMSemi,pasMTrend,pasMAll,pasMTimeProper = seasonalCycle(PasMonthly,PasMonthlyTime)
    pvMSeason,pvMSemi,pvMTrend,pvMAll,pvMTimeProper = seasonalCycle(PVMonthly,PVMonthlyTime)

#continues from last sanity check but these are just the plots
if sanityCheckDecompositionMonthlyDataPlot:
    genericPlot(pasMTimeProper,pasMSeason,"pas monthly avg season")
    genericPlot(pasMTimeProper,pasMAll,"pas monthly avg all")
    genericPlot(pasMTimeProper,pasMSemi,"pas monthly avg semi")
    genericPlot(pasMTimeProper,pasMTrend,"pas monthly avg trend")
    
    genericPlot(pvMTimeProper,pvMSeason,"pv monthly avg season")
    genericPlot(pvMTimeProper,pvMAll,"pv monthly avg all")
    genericPlot(pvMTimeProper,pvMSemi,"pv monthly avg semi")
    genericPlot(pvMTimeProper,pvMTrend,"pv monthly avg trend")

### sanity check-- reconstruction of monthly data from decomposition - plot with original data
if sanityCheckReconstructionMonthlyData:
    #reconstruction
    reconstructedPas = reconstruct(pasMSeason,pasMSemi,pasMTrend,pasMAll)
    reconstructedPV = reconstruct(pvMSeason,pvMSemi,pvMTrend,pvMAll)
    
    """
    genericPlot(pasMTimeProper,reconstructedPas,"pas reconstructed")
    genericPlot(pvMTimeProper,reconstructedPV,"pv reconstructed")
    """
    
    decompOriginPlot(pasMTimeProper,reconstructedPas,PasMonthlyTime,PasMonthly,"pasadena")
    decompOriginPlot(pvMTimeProper,reconstructedPV,PVMonthlyTime,PVMonthly,"palos verdes")
    
### sanity check-- plot reconstruction of monthly data + aspects
if sanityCheckReconstructionAspects:
    seasonalCyclePlot(PasMonthly,PasMonthlyTime,"pasadena")
    seasonalCyclePlot(PVMonthly,PVMonthlyTime,"palos verdes")

### sanity check-- average annual cycle decomposed from scratch
if sanityCheckAnnualDecomposedScratch:
    # PV - 1. shave datasets
    pvCO2_1,pvTime_1 = shave(pvCO2,pvTime)
    pvCO2_2 = pvCO2_1
    # PAS - 1. doMath (interpolates NaN, calcs trendline, removes trendline) FOR PASADENA, remove before 2004
    pasCO2_2,pasTime_2 = removeBefore2004(pasCO2,pasTime)
    pasCO2_2_5 = doMath(np.array(pasCO2_2))
    # PV - 2. calculate and remove trendline for palos verdes, ALSO REMOVE NAN VALUE FROM PV
    pvCO2_2_5,pvTime_2 = removeNaN(list(pvCO2_2),list(pvTime_1))
    pvCO2_2_5 = np.array(pvCO2_2_5)
    pvTime_2 = pandas.core.series.Series(pvTime_2)
    m_pv,b_pv = bestFit(pvCO2_2_5)
    pvCO2_3 = removeTrendline(m_pv,b_pv,pvCO2_2_5)
    # 2. make datetime stamps for one year for pasadena
    yearDT_PAS = makeDTYear(2001,"pandas")
    # 2. make datetime stamps for one week for palos verdes
    yearDT_PV = makeDTYear_m(2001,"pandas")
    """
    # 3. average data into one year
    pasCO2_3 = averageYears_day(pandas.DataFrame({'date':pasTime_2,'co2':pasCO2_2_5},columns=['date','co2']))
    pvCO2_4 = averageYears_month(pandas.DataFrame({'date':pvTime_1,'co2':pvCO2_3},columns=['date','co2']))
    # 4. seasonal cycle the avg data
    # return (season,semi,trend,all_list,convertedTimeTemp)
    result = removeNaN(pvCO2_4,list(yearDT_PV))
    pvCO2_4 = result[0]
    yearDT_PV = pandas.core.series.Series(result[1])
    pasSeason_4,pasSemi_4,pasTrend_4,pasAll_4,pasTime_4 = seasonalCycle(pasCO2_3,yearDT_PAS)
    pvSeason_4,pvSemi_4,pvTrend_4,pvAll_4,pvTime_4 = seasonalCycle(pvCO2_4,yearDT_PV)
    # 5. plot everything
    """
    """
    genericPlot(pasTime_4,pasSeason_4,"pas yearly season")
    genericPlot(pasTime_4,pasSemi_4,"pas yearly semi")
    genericPlot(pasTime_4,pasTrend_4,"pas yearly trend")
    genericPlot(pasTime_4,pasAll_4,"pas yearly all")
    
    genericPlot(pvTime_4,pvSeason_4,"pv yearly season")
    genericPlot(pvTime_4,pvSemi_4,"pv yearly semi")
    genericPlot(pvTime_4,pvTrend_4,"pv yearly trend")
    genericPlot(pvTime_4,pvAll_4,"pv yearly all")
    
    plt.figure()
    plt.plot(yearDT_PV,pvCO2_4,'blue')
    plt.plot(pvTime_4,pvAll_4,'green')
    plt.legend(handles=[pc.Patch(color='green',label='decomposed'),pc.Patch(color='blue',label='data')])
    plt.title("palos verdes, decomposed vs. actual")
    plt.show()
    
    plt.figure()
    plt.plot(yearDT_PAS,pasCO2_3,'blue')
    plt.plot(pasTime_4,pasAll_4,'green')
    plt.legend(handles=[pc.Patch(color='green',label='decomposed'),pc.Patch(color='blue',label='data')])
    plt.title("pasadena, decomposed vs. actual")
    plt.show()
    """

### 2012 drought split
if droughtSplit:
    #pvCO2_3, pvTime_1
    #pasCO2_2_5,pasTime_2
    """
    result = removeNaN(list(pvCO2_3),list(pvTime_1))
    pvCO2_3 = np.array(result[0])
    pvTime_1 = pandas.core.series.Series(result[1])
    """
    
    pvSplit2013 = splitBeforeYear(pvCO2_3,pvTime_2,2013)
    pvSplit2013After = pvSplit2013[0]
    pvSplit2013Before = pvSplit2013[1]
    
    #for general comparison purposes
    pvCO2_4,pvTime_3 = removeDataPoint(pvCO2_3,pvTime_2,pandas.Timestamp('2014-08-30 00:00:00'))
    
    #take 8-30-2014 data out of pvSplit2013[1]
    pvDataAfter_fixed,pvTimeAfter_fixed = removeDataPoint(pvSplit2013After[0],pvSplit2013After[1],pandas.Timestamp('2014-08-30 00:00:00'))
    
    pasSplit2013 = splitBeforeYear(pasCO2_2_5,pasTime_2,2013)
    pasSplit2013After = pasSplit2013[0]
    pasSplit2013Before = pasSplit2013[1]

### plotting that accompanies droughtSplit code
if droughtSplitPlot:
    
    seasonalCyclePlot_modified(pvSplit2013Before[0],pvSplit2013Before[1],"PV through 2012")
    seasonalCyclePlot_modified(pvDataAfter_fixed,pvTimeAfter_fixed,"PV after 2012")
    
    seasonalCyclePlot_modified(pasSplit2013Before[0],pasSplit2013Before[1],"Pas through 2012")
    seasonalCyclePlot_modified(pasSplit2013After[0],pasSplit2013After[1],"Pas after 2012")

### average into annual cycle
if droughtSplitAnnual: 
    pvSplit2013 = splitBeforeYear(pvCO2_3,pvTime_2,2013)
    pvSplit2013After = pvSplit2013[0]
    pvSplit2013Before = pvSplit2013[1]
    pvB2013_season,pvB2013_semi,pvB2013_trend,pvB2013_all,pvB2013_time = seasonalCycle(pvSplit2013Before[0],pvSplit2013Before[1]) #dset,time
    pvA2013_season,pvA2013_semi,pvA2013_trend,pvA2013_all,pvA2013_time = seasonalCycle(pvDataAfter_fixed,pvTimeAfter_fixed)
    
    pasSplit2013 = splitBeforeYear(pasCO2_2_5,pasTime_2,2013)
    pasSplit2013After = pasSplit2013[0]
    pasSplit2013Before = pasSplit2013[1]
    pasB2013_season,pasB2013_semi,pasB2013_trend,pasB2013_all,pasB2013_time = seasonalCycle(pasSplit2013Before[0],pasSplit2013Before[1])
    pasA2013_season,pasA2013_semi,pasA2013_trend,pasA2013_all,pasA2013_time = seasonalCycle(pasSplit2013After[0],pasSplit2013After[1])
    
    #12-29-17 note: what I'm working with (one year of pasadena and palos verdes seasonal and semiannual cycles) has to come from the above mentioned variables in this if clause
    
    #pasBefore_seasonAvg = averageYears_month(pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_season},columns=['date','co2']))
    pasBefore_seasonListByMonth = [annualPlotMonthPoint(pasB2013_season,pasB2013_time,x) for x in range(1,13)]
    pasBefore_seasonAvg = [np.mean(pasBefore_seasonListByMonth[x]) for x in range(0,12)]
    pasBefore_seasonStd = [np.std(pasBefore_seasonListByMonth[x]) for x in range(0,12)]
    
    #pasBefore_semiAvg = averageYears_month(pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_semi},columns=['date','co2']))
    pasBefore_semiListByMonth = [annualPlotMonthPoint(pasB2013_semi,pasB2013_time,x) for x in range(1,13)]
    pasBefore_semiAvg = [np.mean(pasBefore_semiListByMonth[x]) for x in range(0,12)]
    pasBefore_semiStd = [np.std(pasBefore_semiListByMonth[x]) for x in range(0,12)]
    
    #pasBefore_trendAvg = averageYears_month(pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_trend},columns=['date','co2']))
    pasBefore_trendListByMonth = [annualPlotMonthPoint(pasB2013_trend,pasB2013_time,x) for x in range(1,13)]
    pasBefore_trendAvg = [np.mean(pasBefore_trendListByMonth[x]) for x in range(0,12)]
    pasBefore_trendStd = [np.std(pasBefore_trendListByMonth[x]) for x in range(0,12)]
    
    #pasBefore_allAvg = averageYears_month(pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_all},columns=['date','co2']))
    pasBefore_allListByMonth = [annualPlotMonthPoint(pasB2013_all,pasB2013_time,x) for x in range(1,13)]
    pasBefore_allAvg = [np.mean(pasBefore_allListByMonth[x]) for x in range(0,12)]
    pasBefore_allStd = [np.std(pasBefore_allListByMonth[x]) for x in range(0,12)]
    
    #pasAfter_seasonAvg = averageYears_month(pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_season},columns=['date','co2']))
    pasAfter_seasonListByMonth = [annualPlotMonthPoint(pasA2013_season,pasA2013_time,x) for x in range(1,13)]
    pasAfter_seasonAvg = [np.mean(pasAfter_seasonListByMonth[x]) for x in range(0,12)]
    pasAfter_seasonStd = [np.std(pasAfter_seasonListByMonth[x]) for x in range(0,12)]
    
    #pasAfter_semiAvg = averageYears_month(pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_semi},columns=['date','co2']))
    pasAfter_semiListByMonth = [annualPlotMonthPoint(pasA2013_semi,pasA2013_time,x) for x in range(1,13)]
    pasAfter_semiAvg = [np.mean(pasAfter_semiListByMonth[x]) for x in range(0,12)]
    pasAfter_semiStd = [np.std(pasAfter_semiListByMonth[x]) for x in range(0,12)]
    
    #pasAfter_trendAvg = averageYears_month(pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_trend},columns=['date','co2']))
    pasAfter_trendListByMonth = [annualPlotMonthPoint(pasA2013_trend,pasA2013_time,x) for x in range(1,13)]
    pasAfter_trendAvg = [np.mean(pasAfter_trendListByMonth[x]) for x in range(0,12)]
    pasAfter_trendStd = [np.std(pasAfter_trendListByMonth[x]) for x in range(0,12)]
    
    #pasAfter_allAvg = averageYears_month(pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_all},columns=['date','co2']))
    pasAfter_allListByMonth = [annualPlotMonthPoint(pasA2013_all,pasA2013_time,x) for x in range(1,13)]
    pasAfter_allAvg = [np.mean(pasAfter_allListByMonth[x]) for x in range(0,12)]
    pasAfter_allStd = [np.std(pasAfter_allListByMonth[x]) for x in range(0,12)]

    #PV DATA WORK
    """
    pvBefore_seasonAvg = averageYears_month(pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_season},columns=['date','co2']))
    pvBefore_semiAvg = averageYears_month(pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_semi},columns=['date','co2']))
    pvBefore_trendAvg = averageYears_month(pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_trend},columns=['date','co2']))
    pvBefore_allAvg = averageYears_month(pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_all},columns=['date','co2']))
    pvAfter_seasonAvg = averageYears_month(pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_season},columns=['date','co2']))
    pvAfter_semiAvg = averageYears_month(pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_semi},columns=['date','co2']))
    pvAfter_trendAvg = averageYears_month(pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_trend},columns=['date','co2']))
    pvAfter_allAvg = averageYears_month(pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_all},columns=['date','co2']))
    """
    
    pvBefore_seasonListByMonth = [annualPlotMonthPoint(pvB2013_season,pvB2013_time,x) for x in range(1,13)]
    pvBefore_seasonAvg = [np.mean(pvBefore_seasonListByMonth[x]) for x in range(0,12)]
    pvBefore_seasonStd = [np.std(pvBefore_seasonListByMonth[x]) for x in range(0,12)]
    
    pvBefore_semiListByMonth = [annualPlotMonthPoint(pvB2013_semi,pvB2013_time,x) for x in range(1,13)]
    pvBefore_semiAvg = [np.mean(pvBefore_semiListByMonth[x]) for x in range(0,12)]
    pvBefore_semiStd = [np.std(pvBefore_semiListByMonth[x]) for x in range(0,12)]
    
    pvBefore_trendListByMonth = [annualPlotMonthPoint(pvB2013_trend,pvB2013_time,x) for x in range(1,13)]
    pvBefore_trendAvg = [np.mean(pvBefore_trendListByMonth[x]) for x in range(0,12)]
    pvBefore_trendStd = [np.std(pvBefore_trendListByMonth[x]) for x in range(0,12)]
    
    pvBefore_allListByMonth = [annualPlotMonthPoint(pvB2013_all,pvB2013_time,x) for x in range(1,13)]
    pvBefore_allAvg = [np.mean(pvBefore_allListByMonth[x]) for x in range(0,12)]
    pvBefore_allStd = [np.std(pvBefore_allListByMonth[x]) for x in range(0,12)]
    
    pvAfter_seasonListByMonth = [annualPlotMonthPoint(pvA2013_season,pvA2013_time,x) for x in range(1,13)]
    pvAfter_seasonAvg = [np.mean(pvAfter_seasonListByMonth[x]) for x in range(0,12)]
    pvAfter_seasonStd = [np.std(pvAfter_seasonListByMonth[x]) for x in range(0,12)]
    
    pvAfter_semiListByMonth = [annualPlotMonthPoint(pvA2013_semi,pvA2013_time,x) for x in range(1,13)]
    pvAfter_semiAvg = [np.mean(pvAfter_semiListByMonth[x]) for x in range(0,12)]
    pvAfter_semiStd = [np.std(pvAfter_semiListByMonth[x]) for x in range(0,12)]
    
    pvAfter_trendListByMonth = [annualPlotMonthPoint(pvA2013_trend,pvA2013_time,x) for x in range(1,13)]
    pvAfter_trendAvg = [np.mean(pvAfter_trendListByMonth[x]) for x in range(0,12)]
    pvAfter_trendStd = [np.std(pvAfter_trendListByMonth[x]) for x in range(0,12)]
    
    pvAfter_allListByMonth = [annualPlotMonthPoint(pvA2013_all,pvA2013_time,x) for x in range(1,13)]
    pvAfter_allAvg = [np.mean(pvAfter_allListByMonth[x]) for x in range(0,12)]
    pvAfter_allStd = [np.std(pvAfter_allListByMonth[x]) for x in range(0,12)]
    
    pasBeforeReconstructedAvg = reconstruct(pasBefore_allAvg,pasBefore_seasonAvg,pasBefore_semiAvg,pasBefore_trendAvg)
    
    pvBeforeReconstructedAvg = reconstruct(pvBefore_allAvg,pvBefore_seasonAvg,pvBefore_semiAvg,pvBefore_trendAvg)
    
    pasAfterReconstructedAvg = reconstruct(pasAfter_allAvg,pasAfter_seasonAvg,pasAfter_semiAvg,pasAfter_trendAvg)
    
    pvAfterReconstructedAvg = reconstruct(pvAfter_allAvg,pvAfter_seasonAvg,pvAfter_semiAvg,pvAfter_trendAvg)
    
    """
    pasBeforeDataAvg = averageYears_month(pandas.DataFrame({'date':pasSplit2013Before[1],'co2':pasSplit2013Before[0]},columns=['date','co2']))
    pvBeforeDataAvg = averageYears_month(pandas.DataFrame({'date':pvSplit2013Before[1],'co2':pvSplit2013Before[0]},columns=['date','co2']))
    pasAfterDataAvg = averageYears_month(pandas.DataFrame({'date':pasSplit2013After[1],'co2':pasSplit2013After[0]},columns=['date','co2']))
    pvAfterDataAvg = averageYears_month(pandas.DataFrame({'date':pvSplit2013After[1],'co2':pvSplit2013After[0]},columns=['date','co2']))
    """
    
    #NEW AND IMPROVED METHOD - 9-18-17 - REPLACES ABOVE COMMENTED CHUNK
    l1 = []
    l2 = []
    l3 = []
    l4 = []
    for x in range (1,13):
        l1.append(annualPlotMonthPoint(pasSplit2013Before[0],pasSplit2013Before[1],x))
        l2.append(annualPlotMonthPoint(pasSplit2013After[0],pasSplit2013After[1],x))
        l3.append(annualPlotMonthPoint(pvSplit2013Before[0],pvSplit2013Before[1],x))
        l4.append(annualPlotMonthPoint(pvDataAfter_fixed,pvTimeAfter_fixed,x))
    pasBeforeDataAvg = [l1[x][0] for x in range (0,len(l1))]
    pasBeforeDataList = [l1[x][1] for x in range (0,len(l1))]
    pvBeforeDataAvg = [l3[x][0] for x in range (0,len(l3))]
    pvBeforeDataList = [l3[x][1] for x in range (0,len(l3))]
    pasAfterDataAvg = [l2[x][0] for x in range (0,len(l2))]
    pasAfterDataList = [l2[x][1] for x in range (0,len(l2))]
    pvAfterDataAvg = [l4[x][0] for x in range (0,len(l4))]
    pvAfterDataList = [l4[x][1] for x in range (0,len(l4))]
    
### plotting that accompanies droughtSplitAnnual
if droughtSplitAnnualPlot:   
    decompOriginPlot(range(1,13),pasBeforeReconstructedAvg,range(1,13),pasBeforeDataAvg,"pas reconstructed vs data, through 2012, avg")
    decompOriginPlot(range(1,13),pvBeforeReconstructedAvg,range(1,13),pvBeforeDataAvg,"pv reconstructed vs data, through 2012, avg")
    decompOriginPlot(range(1,13),pasAfterReconstructedAvg,range(1,13),pasAfterDataAvg,"pas reconstructed vs data, after 2012, avg")
    decompOriginPlot(range(1,13),pvAfterReconstructedAvg,range(1,13),pvAfterDataAvg,"pv reconstructed vs data, after 2012, avg")
    #COMPARISONS B/T BEFORE AND DURING DROUGHT
    plt.figure()
    plt.plot(range(1,13),pasBeforeReconstructedAvg,'--',color='blue')
    plt.plot(range(1,13),pasBeforeDataAvg,'red')
    plt.plot(range(1,13),pasAfterReconstructedAvg,'--',color='green')
    plt.plot(range(1,13),pasAfterDataAvg,'purple')
    plt.legend(handles=[pc.Patch(label='before reconstructed',color='blue'),pc.Patch(label='before actual',color='red'),pc.Patch(label = 'during reconstructed',color='green'),pc.Patch(label='during actual',color='purple')])
    plt.title("Pasadena before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,13,1)))
    plt.show()
    
    plt.figure()
    plt.plot(range(1,13),pvBeforeReconstructedAvg,'--',color='blue')
    plt.plot(range(1,13),pvBeforeDataAvg,'red')
    plt.plot(range(1,13),pvAfterReconstructedAvg,'--',color='green')
    plt.plot(range(1,13),pvAfterDataAvg,'purple')
    plt.legend(handles=[pc.Patch(label='before reconstructed',color='blue'),pc.Patch(label='before actual',color='red'),pc.Patch(label = 'during reconstructed',color='green'),pc.Patch(label='during actual',color='purple')])
    plt.title("Palos Verdes before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,13,1)))
    plt.show()

### droughtSplitAnnualPlot, but with errorbars
if droughtSplitAnnualPlotErrorbars:
    #ERRORBAR WORK
    """
    pvSplit2013Before = 0
    pvSplit2013After = 1
    pasSplit2013Before = 2
    pasSplit2013After = 3
    
    If you calculate the standard errors of the original data’s overall averages (standard deviation of the 4 Jan values (2013-2016) divided by the square root of the number of points (4^2)), then we will have some idea as to how good the final fits are.  """
    errorbars = [] #final error lists
    for x in range (0,4):
        if x == 0:
            dataList = pvBeforeDataList
        elif x == 1:
            dataList = pvAfterDataList
        elif x == 2:
            dataList = pasBeforeDataList
        elif x == 3:
            dataList = pasAfterDataList
        #errorbars.append(np.std(dataList)/math.pow(len(dataList),0.5))
        errorPerYear = []
        for y in range (1,13):
            errorPerYear.append(np.std(dataList[y-1])/math.pow(len(dataList[y-1]),0.5))
        errorbars.append(errorPerYear)
            
    decompOriginPlotErrorbars(range(1,13),pasBeforeReconstructedAvg,range(1,13),pasBeforeDataAvg,"pas reconstructed vs data, through 2012, avg",errorbars[2])
    decompOriginPlotErrorbars(range(1,13),pvBeforeReconstructedAvg,range(1,13),pvBeforeDataAvg,"pv reconstructed vs data, through 2012, avg",errorbars[0])
    decompOriginPlotErrorbars(range(1,13),pasAfterReconstructedAvg,range(1,13),pasAfterDataAvg,"pas reconstructed vs data, after 2012, avg",errorbars[3])
    decompOriginPlotErrorbars(range(1,13),pvAfterReconstructedAvg,range(1,13),pvAfterDataAvg,"pv reconstructed vs data, after 2012, avg",errorbars[1])
    
    #COMPARISONS B/T BEFORE AND DURING DROUGHT 
    plt.figure()
    plt.plot(range(1,13),pasBeforeReconstructedAvg,'--',color='blue')
    plt.errorbar(range(1,13),pasBeforeDataAvg,yerr=errorbars[2],color='blue',fmt='o',capsize=3)
    plt.plot(range(1,13),pasAfterReconstructedAvg,'--',color='red')
    plt.errorbar(range(1,13),pasAfterDataAvg,yerr=errorbars[3],color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Pasadena before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,13,1)))
    plt.show()
    
    plt.figure()
    plt.plot(range(1,13),pvBeforeReconstructedAvg,'--',color='blue')
    plt.errorbar(range(1,13),pvBeforeDataAvg,yerr=errorbars[0],color='blue',fmt='o',capsize=3)
    plt.plot(range(1,13),pvAfterReconstructedAvg,'--',color='red')
    plt.errorbar(range(1,13),pvAfterDataAvg,yerr=errorbars[1],color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Palos Verdes before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,13,1)))
    plt.show()

### extended droughtSplitAnnualPlot plots with month axes 1-18
if droughtSplitAnnualExtended: #repeater()
    errorbars = [] #final error lists
    for x in range (0,4):
        if x == 0:
            dataList = pvBeforeDataList
        elif x == 1:
            dataList = pvAfterDataList
        elif x == 2:
            dataList = pasBeforeDataList
        elif x == 3:
            dataList = pasAfterDataList
        #errorbars.append(np.std(dataList)/math.pow(len(dataList),0.5))
        errorPerYear = []
        for y in range (1,13):
            errorPerYear.append(np.std(dataList[y-1])/math.pow(len(dataList[y-1]),0.5))
        errorbars.append(errorPerYear)
            
    decompOriginPlotErrorbars(range(1,19),repeater(pasBeforeReconstructedAvg,18),range(1,19),repeater(pasBeforeDataAvg,18),"pas reconstructed vs data, through 2012, avg",repeater(errorbars[2],18))
    decompOriginPlotErrorbars(range(1,19),repeater(pvBeforeReconstructedAvg,18),range(1,19),repeater(pvBeforeDataAvg,18),"pv reconstructed vs data, through 2012, avg",repeater(errorbars[0],18))
    decompOriginPlotErrorbars(range(1,19),repeater(pasAfterReconstructedAvg,18),range(1,19),repeater(pasAfterDataAvg,18),"pas reconstructed vs data, after 2012, avg",repeater(errorbars[3],18))
    decompOriginPlotErrorbars(range(1,19),repeater(pvAfterReconstructedAvg,18),range(1,19),repeater(pvAfterDataAvg,18),"pv reconstructed vs data, after 2012, avg",repeater(errorbars[1],18))
    
    #COMPARISONS B/T BEFORE AND DURING DROUGHT
    #9/20/17: COLORS BY SITE
    plt.figure()
    plt.plot(range(1,19),repeater(pasBeforeReconstructedAvg,18),'--',color='blue')
    plt.errorbar(range(1,19),repeater(pasBeforeDataAvg,18),yerr=repeater(errorbars[2],18),color='blue',fmt='o',capsize=3)
    plt.plot(range(1,19),repeater(pasAfterReconstructedAvg,18),'--',color='red')
    plt.errorbar(range(1,19),repeater(pasAfterDataAvg,18),yerr=repeater(errorbars[3],18),color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Pasadena reconstruction (dashed) and average data (point)")
    plt.xlabel("Time (month)")
    plt.ylabel("CO2 concentration (ppm)")
    plt.xticks(list(np.arange(1,19,1)))
    plt.show()
    
    plt.figure()
    plt.plot(range(1,19),repeater(pvBeforeReconstructedAvg,18),'--',color='blue')
    plt.errorbar(range(1,19),repeater(pvBeforeDataAvg,18),yerr=repeater(errorbars[0],18),color='blue',fmt='o',capsize=3)
    plt.plot(range(1,19),repeater(pvAfterReconstructedAvg,18),'--',color='red')
    plt.errorbar(range(1,19),repeater(pvAfterDataAvg,18),yerr=repeater(errorbars[1],18),color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Palos Verdes reconstruction (dashed) and average data (point)")
    plt.xlabel("Time (month)")
    plt.ylabel("CO2 concentration (ppm)")
    plt.xticks(list(np.arange(1,19,1)))
    plt.show()

### truncated droughtSplitAnnualPlot plots with month axes 1-6
if droughtSplitAnnualTruncated: #repeater()
    errorbars = [] #final error lists
    for x in range (0,4):
        if x == 0:
            dataList = pvBeforeDataList
        elif x == 1:
            dataList = pvAfterDataList
        elif x == 2:
            dataList = pasBeforeDataList
        elif x == 3:
            dataList = pasAfterDataList
        #errorbars.append(np.std(dataList)/math.pow(len(dataList),0.5))
        errorPerYear = []
        for y in range (1,13):
            errorPerYear.append(np.std(dataList[y-1])/math.pow(len(dataList[y-1]),0.5))
        errorbars.append(errorPerYear)
            
    decompOriginPlotErrorbars(range(1,7),repeater(pasBeforeReconstructedAvg,6),range(1,7),repeater(pasBeforeDataAvg,6),"pas reconstructed vs data, through 2012, avg",repeater(errorbars[2],6))
    decompOriginPlotErrorbars(range(1,7),repeater(pvBeforeReconstructedAvg,6),range(1,7),repeater(pvBeforeDataAvg,6),"pv reconstructed vs data, through 2012, avg",repeater(errorbars[0],6))
    decompOriginPlotErrorbars(range(1,7),repeater(pasAfterReconstructedAvg,6),range(1,7),repeater(pasAfterDataAvg,6),"pas reconstructed vs data, after 2012, avg",repeater(errorbars[3],6))
    decompOriginPlotErrorbars(range(1,7),repeater(pvAfterReconstructedAvg,6),range(1,7),repeater(pvAfterDataAvg,6),"pv reconstructed vs data, after 2012, avg",repeater(errorbars[1],6))
    
    #COMPARISONS B/T BEFORE AND DURING DROUGHT
    #9/20/17: COLORS BY SITE
    plt.figure()
    plt.plot(range(1,7),repeater(pasBeforeReconstructedAvg,6),'--',color='blue')
    plt.errorbar(range(1,7),repeater(pasBeforeDataAvg,6),yerr=repeater(errorbars[2],6),color='blue',fmt='o',capsize=3)
    plt.plot(range(1,7),repeater(pasAfterReconstructedAvg,6),'--',color='red')
    plt.errorbar(range(1,7),repeater(pasAfterDataAvg,6),yerr=repeater(errorbars[3],6),color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Pasadena before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,7,1)))
    plt.show()
    
    plt.figure()
    plt.plot(range(1,7),repeater(pvBeforeReconstructedAvg,6),'--',color='blue')
    plt.errorbar(range(1,7),repeater(pvBeforeDataAvg,6),yerr=repeater(errorbars[0],6),color='blue',fmt='o',capsize=3)
    plt.plot(range(1,7),repeater(pvAfterReconstructedAvg,6),'--',color='red')
    plt.errorbar(range(1,7),repeater(pvAfterDataAvg,6),yerr=repeater(errorbars[1],6),color='red',fmt='o',capsize=3)
    plt.legend(handles=[pc.Patch(label='before drought',color='blue'),pc.Patch(label='during drought',color='red')])
    plt.title("Palos Verdes before/during drought")
    genericAxes()
    plt.xticks(list(np.arange(1,7,1)))
    plt.show()

### log normal, drought split
if logNormalDroughtSplit:
    plt.figure()
    logNormal(pasSplit2013Before[0],-50,100,1000,"Pasadena Before Drought",'blue','linear','linspace')
    logNormal(pasSplit2013After[0],-50,100,1000,"Pasadena During Drought",'red','linear','linspace')
    plt.legend(handles=[pc.Patch(color = 'blue',label = 'Before Drought'),pc.Patch(color='red',label='During Drought')])
    plt.title("Pasadena, drought split log normal")
    plt.show()
    
    plt.figure()
    logNormal(pvSplit2013Before[0],-50,100,1000,"PV Before Drought",'green','linear','linspace')
    logNormal(pvSplit2013After[0],-50,100,1000,"PV During Drought",'orange','linear','linspace')
    plt.legend(handles=[pc.Patch(color = 'green',label = 'Before Drought'),pc.Patch(color='orange',label='During Drought')])
    plt.title("Palos Verdes, drought split log normal")
    plt.show()
    
### log normal and histograms, drought split
if logNormalHistogramsDroughtSplit:
    plt.figure()
    histogram(pasSplit2013Before[0],100,'blue')
    logNormal(pasSplit2013Before[0],-50,100,1000,"Pasadena Before Drought",'red','linear','linspace')
    plt.title("Pasadena Before Drought log normal")
    plt.show()
    
    plt.figure()
    histogram(pasSplit2013After[0],100,'blue')
    logNormal(pasSplit2013After[0],-50,100,1000,"Pasadena During Drought",'red','linear','linspace')
    plt.title("Pasadena After Drought log normal")
    plt.show()
    
    plt.figure()
    histogram(pvSplit2013Before[0],100,'orange')
    logNormal(pvSplit2013Before[0],-50,100,1000,"PV Before Drought",'green','linear','linspace')
    plt.title("PV Before Drought log normal")
    plt.show()
    
    plt.figure()
    histogram(pvSplit2013After[0],100,'green')
    logNormal(pvSplit2013After[0],-50,100,1000,"PV During Drought",'orange','linear','linspace')
    plt.title("PV After Drought log normal")
    plt.show()
    
### MAUNA LOA CODE ========================

### before/during drought split setup
if ml_droughtSetup:
    mlCO2_detrendlined = doMath(mlCO2)
    mlBeforeDuring = splitBeforeYear(mlCO2_detrendlined,mlDate,2013)
    mlDuringDrought = mlBeforeDuring[0]
    mlBeforeDrought = mlBeforeDuring[1]
    
### log normal drought split
    ### before/during plots together
    if ml_droughtLogNormalTogether:
        plt.figure()
        logNormal(mlBeforeDrought[0],-20,20,1000,"ML Before Drought",'green','linear','linspace')
        logNormal(mlDuringDrought[0],-20,20,1000,"ML During Drought",'orange','linear','linspace')
        plt.legend(handles=[pc.Patch(label='Before drought',color='green'),pc.Patch(label='During drought',color='orange')])
        plt.title("Mauna Loa, Before/During Drought Log Normal")
        plt.show()
        
    ### before/during plots separate with histograms
    if ml_droughtLogNormalHistograms:
        plt.figure()
        histogram(mlBeforeDrought[0],100,'orange')
        logNormal(mlBeforeDrought[0],-20,20,1000,"ML Before Drought",'green','linear','linspace')
        plt.title("Mauna Loa before drought log normal")
        plt.show()
        
        plt.figure()
        histogram(mlDuringDrought[0],100,'orange')
        logNormal(mlDuringDrought[0],-20,20,1000,"ML During Drought",'green','linear','linspace')
        plt.title("Mauna Loa During drought log normal")
        plt.show()

### MAUNA LOA AND PV/PAS=====================
"""
pv: pvCO2_3, pvTime_2
pas: pasCO2_2_5, pasTime_2
ml: mlCO2_detrended, mlDate
"""
if ml_pv_pas_logNormal:
    plt.figure()
    logNormal(pvCO2_3,-40,60,1000,"",'green','linear','linspace')
    logNormal(pasCO2_2_5,-40,60,1000,"",'red','linear','linspace')
    logNormal(mlCO2_detrendlined,-40,60,1000,"",'blue','linear','linspace')
    plt.legend(handles=[pc.Patch(label='Palos Verdes',color='green'),\
                        pc.Patch(label='Pasadena',color='red'),\
                        pc.Patch(label='Mauna Loa',color='blue')])
    plt.title("Pasadena, Mauna Loa, Palos Verdes log normal")
    plt.show()

"""
pv: pvSplit2013Before, pvSplit2013After
pas: pasSplit2013Before, pasSplit2013After
ml: mlBeforeDrought, mlDuringDrought
"""
if ml_pv_pas_logNormalSplit:
    #before drought
    plt.figure()
    logNormal(pvSplit2013Before[0],-40,60,1000,"",'green','linear','linspace')
    logNormal(pasSplit2013Before[0],-40,60,1000,"",'red','linear','linspace')
    logNormal(mlBeforeDrought[0],-40,60,1000,"",'blue','linear','linspace')
    plt.legend(handles=[pc.Patch(label='Palos Verdes',color='green'),\
                        pc.Patch(label='Pasadena',color='red'),\
                        pc.Patch(label='Mauna Loa',color='blue')])
    plt.title("Pasadena, Mauna Loa, Palos Verdes log normal - Through 2012")
    plt.show()
    
    #after drought
    plt.figure()
    logNormal(pvSplit2013After[0],-40,60,1000,"",'green','linear','linspace')
    logNormal(pasSplit2013After[0],-40,60,1000,"",'red','linear','linspace')
    logNormal(mlDuringDrought[0],-40,60,1000,"",'blue','linear','linspace')
    plt.legend(handles=[pc.Patch(label='Palos Verdes',color='green'),\
                        pc.Patch(label='Pasadena',color='red'),\
                        pc.Patch(label='Mauna Loa',color='blue')])
    plt.title("Pasadena, Mauna Loa, Palos Verdes log normal - After 2012")
    plt.show()

### last fifteen years data manipulation
if lastFifteenYears_dataManipulation:
    #ml_lastFifteenYears = splitBeforeYear(mlCO2,mlDate,mlDate.iloc[mlDate.size-1].year-15)[0]
    pas_lastFifteenYears = splitBeforeYear(pasCO2,pasTime,2004)[0] #used to be pasTime.iloc[pasTime.size-1].year-15, but that goes past 2004
    ljo_lastFifteenYears = splitBeforeYear(laJollaData,laJollaData["Sample Date.1"],2004)[0][0]
    ljo_lastFifteenYears = pandas.concat([laJollaData[ljo_lastFifteenYears.index[0]-1:ljo_lastFifteenYears.index[0]],ljo_lastFifteenYears]) #add last 2003 entry onto dataset for context when using findAdjacentDates to interpolate
    
    #drop index ml_last and pas_last USING GENERATORS: https://wiki.python.org/moin/Generators
    """ml_lastFifteenYearsGenerator = (ml_lastFifteenYears[x].reset_index(drop=True) for x in range (0,2))
    ml_droppedIndex = (ml_lastFifteenYearsGenerator.next(),ml_lastFifteenYearsGenerator.next())"""
    
    pas_lastFifteenYearsGenerator = (pas_lastFifteenYears[x].reset_index(drop=True) for x in range (0,2))
    pas_droppedIndex = (pas_lastFifteenYearsGenerator.next(),pas_lastFifteenYearsGenerator.next())
    
    ljo_droppedIndex = ljo_lastFifteenYears.reset_index(drop=True)
    
    pv_lastFifteenYears = (pvCO2,pvTime)

### last fifteen years histograms
if lastFifteenYears_histograms:
    
    plt.figure()
    #histogram(removeNaN(list(ml_lastFifteenYears[0])),200,"blue")
    #logNormal(removeNaN(list(ml_lastFifteenYears[0])),350,500,1000,"",'red','linear','linspace')
    plt.xlim([350,500])
    plt.title("Mauna Loa")
    plt.show()
    
    plt.figure()
    histogram(removeNaN(list(pas_lastFifteenYears[0])),200,"blue")
    #logNormal(removeNaN(list(pas_lastFifteenYears[0])),350,500,1000,"",'red','linear','linspace')
    plt.xlim([350,500])
    plt.title("Pasadena")
    plt.show()
    
    plt.figure()
    histogram(removeNaN(list(pv_lastFifteenYears[0])),200,"blue")
    #logNormal(removeNaN(list(pv_lastFifteenYears[0])),350,500,1000,"",'red','linear','linspace')
    plt.xlim([350,500])
    plt.title("Palos Verdes")
    plt.show()
    
### process Pasadena CO and O3 data, do lots of work
if CO_O3_correlation:
    #import CO and O3 data, conjoin
    """
    CO_data = pandas.read_csv("C:\Users\Kyle W\Desktop\CO_PICKDATA_2017-10-22.csv")
    O3_data = pandas.read_csv("C:\Users\Kyle W\Desktop\OZONE_PICKDATA_2017-10-22.csv")
    """
    CO_data_master = pandas.DataFrame()
    O3_data_master = pandas.DataFrame()
    
    ozone_list = os.walk("C:\Users\Kyle W\Desktop\Carbon Monoxide and Ozone Aqmis2 Data\Ozone").next()[2]
    monoxide_list = os.walk("C:\Users\Kyle W\Desktop\Carbon Monoxide and Ozone Aqmis2 Data\Carbon Monoxide").next()[2]
    
    for x in range(0,14):
        CO_data = pandas.read_csv("C:\Users\Kyle W\Desktop\Carbon Monoxide and Ozone Aqmis2 Data\Carbon Monoxide\\" + monoxide_list[x],skipfooter = 11,engine='python')
        CO_data_master = pandas.concat([CO_data_master,CO_data],ignore_index=True)
        
        O3_data = pandas.read_csv("C:\Users\Kyle W\Desktop\Carbon Monoxide and Ozone Aqmis2 Data\Ozone\\" + ozone_list[x],skipfooter = 11,engine='python')
        O3_data_master = pandas.concat([O3_data_master,O3_data],ignore_index=True)
        
    #remove la jolla background from pasadena data
    """VARIABLES TO WORK WITH:
        ljo_droppedIndex
        pas_droppedIndex
        
        CONCENTRATION VARIABLES:
            ljo_droppedIndex["Cx"]
            pas_droppedIndex[0]
        TIME VARIABLES:
            ljo_droppedIndex["Sample Date"]
            pas_droppedIndex[1]"""
    #convert La Jolla Date/Time to decimal
    laJollaTimestampDates = ljo_droppedIndex["Sample Date"].copy()
    laJollaDecimalDates = ljo_droppedIndex["Sample Date"].copy()
    pasRemovedLJO = pas_droppedIndex[0].copy() #pas, but without the la jolla background
    for x in range (0,len(ljo_droppedIndex["Sample Time"])):
        year = int(ljo_droppedIndex["Sample Date"][x].split('/')[2])
        month = int(ljo_droppedIndex["Sample Date"][x].split('/')[0])
        day = int(ljo_droppedIndex["Sample Date"][x].split('/')[1])
        hour = int(ljo_droppedIndex["Sample Time"][x][1:].split(':')[0])
        minute = int(ljo_droppedIndex["Sample Time"][x][1:].split(':')[1])
        laJollaTimestampDates[x] = dt(year,month,day,hour,minute)
        laJollaDecimalDates[x] = toYearFraction(laJollaTimestampDates[x])

    
    #find and subtract la jolla background from pasadena, LINEARLY INTERPOLATE BETWEEN LJ DATA POINTS
    for x in range (0,len(pasRemovedLJO)):
        adjacentDates = findAdjacentDates(pas_droppedIndex[1][x],laJollaTimestampDates)
        poly = np.poly1d(np.polyfit(list(laJollaDecimalDates[adjacentDates.index[0]:adjacentDates.index[0]+2]),list(ljo_droppedIndex["Cx"][adjacentDates.index[0]:adjacentDates.index[0]+2]),1))
        pasRemovedLJO[x] -= poly(toYearFraction(pas_droppedIndex[1][x]))
    
    
    #find >2/3 sd's of pasadena distribution
    pasRemovedMean = np.mean(pasRemovedLJO)
    pasRemovedSTD = np.std(pasRemovedLJO)
    pasUpperBound = pasRemovedMean + 2 * pasRemovedSTD #2 sd's from the mean in either direction
    #pasLowerBound = pasRemovedMean - 2 * pasRemovedSTD ||| this is negative, we don't need it
    
    unusualDataIndices = [x for x in range(0,len(pasRemovedLJO)) if pasRemovedLJO[x] > pasUpperBound]
    unusualData = [pasRemovedLJO[x] for x in unusualDataIndices]
    
    #find average of selected pasadena days
    unusualDataAvg = np.mean(unusualData)
    
    #calculate extra CO2 (value - average)
    extraCO2 = [(unusualData[x] - unusualDataAvg) for x in range(0,len(unusualData))]
    extraTime = [pas_droppedIndex[1][x] for x in unusualDataIndices]
    
    #remove seasonal cycle from extra CO2
    """return (season,semi,trend,all_list,convertedTimeTemp)"""
    xSeason,xSemi,xTrend,xAll,xTime = seasonalCycle(extraCO2,extraTime)
    
    extraCO2_1 = [(extraCO2[x]-xSeason[x]) for x in range(0,len(extraCO2))]
    
    """VARIABLES TO WORK WITH
        CO_data_master
        O3_data_master"""
    #fuse CO and O3 measurements taken on the same day
    #CO_fused = pandas.DataFrame(columns=["date","value"])
    #O3_fused = pandas.DataFrame(columns=["date","value"])
    
    #sets have only distinct objects
    CO_dates = list(set([(CO_data_master["date"][x]) for x in range(0,len(CO_data_master))]))
    O3_dates = list(set([(O3_data_master["date"][x]) for x in range(0,len(O3_data_master))]))
    CO_dates.sort()
    O3_dates.sort()
    CO_dates_pandas = [pandas.Timestamp(x) for x in CO_dates]
    O3_dates_pandas = [pandas.Timestamp(x) for x in O3_dates]
    CO_averaged = []
    O3_averaged = []
    for x in range(0,len(CO_dates)):
        CO_averaged.append(np.mean([CO_data_master.iloc[y]["value"] for y in CO_data_master[CO_data_master["date"]==CO_dates[x]].index]))
    for x in range (0,len(O3_dates)):
        O3_averaged.append(np.mean([O3_data_master.iloc[y]["value"] for y in O3_data_master[O3_data_master["date"]==O3_dates[x]].index]))
    
    #this stuff is overkill, leave it alone
    #relevant_indices_CO = [i for i, x in enumerate(CO_dates_pandas) if x == "whatever"]
    #relevant_indices_O3 = [i for i, x in enumerate(O3_dates_pandas) if x == "whatever"]
    
    #indices in extraTime whose corresponding times bear no counterpart in the other time list, with respect to CO and O3
    irrelevant_indices_CO = []
    irrelevant_indices_O3 = []
    
    
    #with respect to carbon monoxide
    for x in range (0,len(extraTime)):
        try:
            CO_dates_pandas.index(extraTime[x])
        except ValueError:
            irrelevant_indices_CO.append(extraTime.index(extraTime[x]))
    
    #just iterate over fresh list to create new sanitized copy
    extraTime_wrt_CO = []
    extraCO2_wrt_CO = []
    for x in range (0,len(extraTime)):
        if x not in irrelevant_indices_CO:
            extraTime_wrt_CO.append(extraTime[x])
            extraCO2_wrt_CO.append(extraCO2_1[x])
        
    #TODO - THE FOLLOWING CODE BLOCK CAN BE MADE INTO A GENERATOR
    CO_averaged_desirable = []
    #CO_time_desirable = [] should theoretically be THE EXACT SAME as extraTime_wrt_CO
    for x in range (0,len(extraTime_wrt_CO)):
        CO_averaged_desirable.append(CO_averaged[CO_dates_pandas.index(extraTime_wrt_CO[x])])

    
    #with respect to ozone
    for x in range (0,len(extraTime)):
        try:
            O3_dates_pandas.index(extraTime[x])
        except ValueError:
            irrelevant_indices_O3.append(extraTime.index(extraTime[x]))
    extraTime_wrt_O3 = []
    extraCO2_wrt_O3 = []
    for x in range (0,len(extraTime)):
        if x not in irrelevant_indices_O3:
            extraTime_wrt_O3.append(extraTime[x])
            extraCO2_wrt_O3.append(extraCO2_1[x])
    O3_averaged_desirable = []
    for x in range (0,len(extraTime_wrt_O3)):
        O3_averaged_desirable.append(O3_averaged[O3_dates_pandas.index(extraTime_wrt_O3[x])])
    
    ### calculate correlation
    #between CO2 and CO
    CO_r = pearsonr(extraCO2_wrt_CO,CO_averaged_desirable)
    
    #between CO2 and O3
    O3_r = pearsonr(extraCO2_wrt_O3,O3_averaged_desirable)
    
if plot_CO_O3_extraCO2:
    #plot CO
    plt.figure()
    plt.scatter(extraTime_wrt_CO,extraCO2_wrt_CO,color="red")
    plt.scatter(extraTime_wrt_CO,CO_averaged_desirable,color="blue")
    plt.xlabel("time")
    plt.ylabel("concentration (ppm)")
    plt.legend(handles=[pc.Patch(color="red",label="extra CO2"),pc.Patch(color="blue",label="CO")])
    plt.show()
    #plot O3
    plt.figure()
    plt.scatter(extraTime_wrt_O3,extraCO2_wrt_O3,color="red")
    plt.scatter(extraTime_wrt_O3,O3_averaged_desirable,color="blue")
    plt.xlabel("time")
    plt.ylabel("concentration (ppm)")
    plt.legend(handles=[pc.Patch(color="red",label="extra CO2"),pc.Patch(color="blue",label="O3")])
    plt.show()   
    
if final_pasadena_three_step:
    scatterplot = False
    regularplot = False
    plot_season_semi_avg = False
    
    #1. original - pas_droppedIndex
    #2. linearly detrended - pas_removed_trendline
    #3. linearly detrended and then deseasonalized - pas_removed_trendseason - (it was never actually deseasonalized)
    #4. plot one year, averaged
    
    #-------------------------------------------------------------------2004 - 2017
    pas_removed_trendline = pas_droppedIndex[0].copy()
    
    #12/4/17 new method: remove linear trend from pas_droppedIndex
    poly = np.poly1d(np.polyfit(pas_droppedIndex[0],[toYearFraction(x) for x in pas_droppedIndex[1]],1))
    for x in range (0,len(pasRemovedLJO)):
        pas_removed_trendline[x] -= poly(toYearFraction(pas_droppedIndex[1][x]))
    
    pas_removed_trendseason = pas_removed_trendline.copy() #never finished work on finding pas_removed_trendseason
    
    f_season,f_semi,f_trend,f_all_list,f_convertedTimeTemp = seasonalCycle(pas_removed_trendseason,pas_droppedIndex[1])
    
    if scatterplot:
        #scatterplot
        plt.figure()
        plt.scatter(f_convertedTimeTemp,f_season,color='orange')
        plt.scatter(f_convertedTimeTemp,f_semi,color='blue')
        plt.legend(handles=[pc.Patch(label='season',color='orange'),pc.Patch(label='semi',color='blue')])
        plt.title("pasadena decomposed")
        plt.show()
    
    if regularplot:
        #regular plot
        plt.figure()
        plt.plot(f_convertedTimeTemp,f_season,color='orange')
        plt.plot(f_convertedTimeTemp,f_semi,color='blue')
        plt.legend(handles=[pc.Patch(label='season',color='orange'),pc.Patch(label='semi',color='blue')])
        plt.title("pasadena decomposed")
        plt.show()
    
    #use averageYears_day to stick this entire thing into a year
    df_pas_f_season = pandas.DataFrame({'date':f_convertedTimeTemp,'co2':f_season},columns=['date','co2'])
    pas_f_season_year = averageYears_day(df_pas_f_season)
    
    df_pas_f_semi = pandas.DataFrame({'date':f_convertedTimeTemp,'co2':f_semi},columns=['date','co2'])
    pas_f_semi_year = averageYears_day(df_pas_f_semi)
    
    if plot_season_semi_avg:
        plt.figure()
        plt.plot(range(1,366),pas_f_season_year,'red')
        plt.plot(range(1,366),pas_f_semi_year,'blue')
    #-------------------------------------------------------------------
    
    pas_droppedIndex_split = splitBeforeYear(pas_droppedIndex[0],pas_droppedIndex[1],2013)
    
    #-------------------------------------------------------------------2004 - 2012
    pas_firstHalf = pas_droppedIndex_split[0]
    
    
    #-------------------------------------------------------------------

    """
    final notes 12-29-17
    leave this alone, go back to sanityCheckAnnualDecomposedScratch and start from there
    """
"""
def makeDataframe(*args):
    #each arg should be a tuple of (column_name/string,column_data/list/array)
    x = {row.SITE_NAME : row.LOOKUP_TABLE for row in cursor}
"""

#12-29-17
def plotOneYearFit(data):
    #a plotting method specific to the data from oneYearFit
    plt.plot(data)
    plt.xticks([1,32,60,91,121,152,182,213,243,274,304,335],('January','February','March','April','May','June','July','August','September','October','November','December'))
    plt.xlabel("Month")
    plt.ylabel("CO2 concentration (ppm)")

if oneYearFit:
    """VARIABLES TO WORK WITH - B is before, A is after
    pvB2013_season,pvB2013_semi,pvB2013_trend,pvB2013_all,pvB2013_time
    pvA2013_season,pvA2013_semi,pvA2013_trend,pvA2013_all,pvA2013_time
    
    pasB2013_season,pasB2013_semi,pasB2013_trend,pasB2013_all,pasB2013_time
    pasA2013_season,pasA2013_semi,pasA2013_trend,pasA2013_all,pasA2013_time
    """
    #BEFORE---
    df_pas_b_season = pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_season},columns=['date','co2'])
    df_pas_b_semi = pandas.DataFrame({'date':pasB2013_time,'co2':pasB2013_semi},columns=['date','co2'])
    
    df_pv_b_season = pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_season},columns=['date','co2'])
    df_pv_b_semi = pandas.DataFrame({'date':pvB2013_time,'co2':pvB2013_semi},columns=['date','co2'])
    
    df_pas_b_season_year = averageYears_day(df_pas_b_season)
    df_pas_b_semi_year = averageYears_day(df_pas_b_semi)
    df_pv_b_season_year = averageYears_day(df_pv_b_season)
    df_pv_b_semi_year = averageYears_day(df_pv_b_semi)
    
    #plot
    plt.figure()
    plt.plot(df_pas_b_season_year,'red')
    plt.plot(df_pas_b_semi_year,'blue')
    plt.xticks([1,32,60,91,121,152,182,213,243,274,304,335],('January','February','March','April','May','June','July','August','September','October','November','December'),rotation=17)
    plt.xlabel("Month")
    plt.ylabel("CO2 concentration (ppm)")
    plt.title("Pasadena Before Drought, Season and Semiannual")
    plt.legend(handles=[pc.Patch(color='red',label='Season'),pc.Patch(color='blue',label='Semiannual')])

    plt.figure()
    plt.plot(df_pv_b_season_year,'red')
    plt.plot(df_pv_b_semi_year,'blue')
    plt.xticks([1,32,60,91,121,152,182,213,243,274,304,335],('January','February','March','April','May','June','July','August','September','October','November','December'),rotation=17)
    plt.xlabel("Month")
    plt.ylabel("CO2 concentration (ppm)")
    plt.title("Palos Verdes Before Drought, Season and Semiannual")
    plt.legend(handles=[pc.Patch(color='red',label='Season'),pc.Patch(color='blue',label='Semiannual')])
    
    #AFTER---
    df_pas_a_season = pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_season},columns=['date','co2'])
    df_pas_a_semi = pandas.DataFrame({'date':pasA2013_time,'co2':pasA2013_semi},columns=['date','co2'])
    
    df_pv_a_season = pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_season},columns=['date','co2'])
    df_pv_a_semi = pandas.DataFrame({'date':pvA2013_time,'co2':pvA2013_semi},columns=['date','co2'])
    
    df_pas_a_season_year = averageYears_day(df_pas_a_season)
    df_pas_a_semi_year = averageYears_day(df_pas_a_semi)
    df_pv_a_season_year = averageYears_day(df_pv_a_season)
    df_pv_a_semi_year = averageYears_day(df_pv_a_semi)
    
    
    #plot
    plt.figure()
    plt.plot(df_pas_a_season_year,'red')
    plt.plot(df_pas_a_semi_year,'blue')
    plt.xticks([1,32,60,91,121,152,182,213,243,274,304,335],('January','February','March','April','May','June','July','August','September','October','November','December'),rotation=17)
    plt.xlabel("Month")
    plt.ylabel("CO2 concentration (ppm)")
    plt.title("Pasadena During/After Drought, Season and Semiannual")
    plt.legend(handles=[pc.Patch(color='red',label='Season'),pc.Patch(color='blue',label='Semiannual')])

    plt.figure()
    plt.plot(df_pv_a_season_year,'red')
    plt.plot(df_pv_a_semi_year,'blue')
    plt.xticks([1,32,60,91,121,152,182,213,243,274,304,335],('January','February','March','April','May','June','July','August','September','October','November','December'),rotation=17)
    plt.xlabel("Month")
    plt.ylabel("CO2 concentration (ppm)")
    plt.title("Palos Verdes After Drought, Season and Semiannual")
    plt.legend(handles=[pc.Patch(color='red',label='Season'),pc.Patch(color='blue',label='Semiannual')])