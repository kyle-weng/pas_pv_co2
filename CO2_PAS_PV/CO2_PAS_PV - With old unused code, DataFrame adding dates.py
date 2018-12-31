#Pasadena and Palos Verdes CO2 Analysis, Kyle Weng, July-August 2017

#to do: Dr. Newman also sent data on more sites
#also: graph data on top of "all"
#   average year data, show one monthly cycle
#also: data for last drought (2012/13 - 15)


#imports
from matplotlib import pyplot as plt, patches as pc
from scipy.stats import lognorm
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

global pastime,pvtime

#common methods
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
    x = np.arange(len(y))
    m,b = np.polyfit(x,y,1)
    return m,b

def interpolateNaN(y):  #checks pandas excel object for NaN values
                        #   and linearly interpolates missing values
                        #   *so long as the NaN values are not at the front or back
                        #y = data
	zz = y.copy()
	xz = np.where(~np.isfinite(y))
	for z in range (0,len(xz)):
		zz[xz[0][z]] = (y[xz[0][z]-1] + y[xz[0][z]+1]) / 2
	return zz

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
        for x in range (0,len(distro)):
            #distro[x] = distro[x] + 40
            #adding 40 really does not matter for the overall distribution, it literally just shifts it right
            print ""
    else:
        print "That didn't work."
        return

    shape,loc,scale = lognorm.fit(dset)
    
    
    pdf = lognorm.pdf(distro, shape, loc, scale)
    
    plt.plot(distro, pdf,color=e) #formerly ax.plot
    plt.title(d + " PDF with data")

def split(dset):
    index = int(len(dset) / 2)
    return (dset[:index],dset[index:])

def histogram(dset,bins,color_in):
    #dset should be self explanatory
    #bins is the number of bins you want the data to be sorted into
    #color_in is the color of the histogram
    plt.hist(dset,bins,normed=1,color=color_in) #formerly ax.hist

def pc_regress(x_in,y_in): #Fai's IDL code but in Python, not IDL
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
    
    for iterator in range (0,6):
        if results.bse[iterator] != float(0):
            stat.append(math.fabs(a.coef_[iterator]/results.bse[iterator]))
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

def genericAxes():
    plt.ylabel('Time (year)')
    plt.xlabel('CO2 concentration (ppm)')
    
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
    global toBeInserted,df_modified,result,index,repeat,toBeInsertedList,x,z
    
    df_modified = pandas.DataFrame()
    #df_modified = df_in.copy()
    
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
    #now we rebuild the dataframe from the ground up
    """
    #concat method
    z = 0
    for x in range (0,len(toBeInsertedList)+df_in.shape[0]):
        if z < len(toBeInsertedList):
            if x == result[1][z][0]: #if x matches an index to be reinserted
                ab = toBeInsertedList[z] #dataframe chunk that should be inserted at that index
                df_modified = pandas.concat([df_modified.iloc[:x+1],ab,df_modified.iloc[x+1:]])
                df_modified.sort_values('date')
                df_modified.reset_index(drop=True)
                x += result[1][z][1]
                z += 1
    """
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
        """
        else:
            ab = df_in.iloc[x:x+1]
        #df_modified = df_modified.append(ab)
        df_modified = pandas.concat([df_modified[:x+1],toBeInserted,df_modified[x+1:]])
        """
    return df_modified
    
    
    """
    global toBeInserted,df_modified
    
    df_modified = df_in.copy()
    while (not checkConsecutiveDates(df_modified)):
        x += 1
        if (df_modified.iloc[x]['date'] + pandas.Timedelta(days=1)) != (df_modified.iloc[x+1]['date']):
            #insert a new row y number of times as necessary
            y = (df_modified.iloc[x+1]['date'] - df_modified.iloc[x]['date']).days-1
            
            #toBeInserted is a fragment of a dataframe that
            #   will be concatenated with the main one
            toBeInserted = pandas.DataFrame()
            for z in range (1,y+1):
                toBeInserted = pandas.concat([toBeInserted,pandas.DataFrame({'date':[df_modified.iloc[x]['date'] + pandas.Timedelta(days=z)],'co2':None})])
            df_modified = pandas.concat([df_modified[:x+1],toBeInserted,df_modified[x+1:]])
            df_modified = df_modified.sort_values('date')
    return df_modified
    """

def averageYears(dset,time): #returns one year of averaged year data
                            #averageYears(pasAll,pasTimeProper)
    global df,uniqueList,comparer,x,pseudoDT
    data = {'date':time,'co2':dset}
    df = pandas.DataFrame(data,columns=['date','co2'])
   
    #insert empty rows to make date column consecutive
    consecutiveDates = False
    return df
    """
    sys.exit(0)
    
    #is there a leap year in the data?  if so, add 2/29 as a data point
    uniqueList = []
    comparer = dt.now()
    for x in range (0,len(time)):
        if (comparer.year != time[x].year):
            comparer = time[x]
            uniqueList.append(time[x].year)
    leap = False #is there a leap year or not
    for x in range (0,len(uniqueList)):
        if uniqueList[x] % 4 == 0:
            leap = True
            break
    pseudoDT = []
    for month in range (1,13):
        if month == 2 and leap == True:
            temp_year = 2000 #there is literally no significance to these years
        else:
            temp_year = 2001
        for day in range (1,monthrange(temp_year,month)[1]+1):
            pseudoDT.append(str(month)+"/"+str(day))
    """
    
        
    
    #interpolate, then average
    
    #first: interpolate, fill in the blanks for each missing day
    
        
    
    #unique_frame = pandas.DataFrame()

local = True
site = "Palos Verdes and Pasadena"

#site specific code begins here
if site == "Palos Verdes and Pasadena":
    if local == True:
        data = pandas.ExcelFile("C:\Users\Kyle W\Desktop\Pas PV flasks CO2 13C 18O summ for Kyle.xlsx") #local
    else:
        data = pandas.ExcelFile("/home/kaw/Pasadena_PalosVerdes_CO2/Pas PV flasks CO2 13C 18O summ for Kyle.xlsx") #server

    pasadena_flasks = data.parse("Pasadena_flasks")
    pasCO2 = pasadena_flasks["CO2 (ppm)"]
    pasTime = pasadena_flasks["date"]
    
    palosverdes_flasks = data.parse("PV_flasks_aves")
    pvCO2 = palosverdes_flasks["[CO2]B"]
    pvTime = palosverdes_flasks["Date"]

    def removeBefore2004(dset,time):
        for x in range (0,len(time)):
            if str(time[x])[:4] == '2004':
                y = x
                break
        return (dset[y:],time[y:])
    
    def seasonalLegend():
        plt.legend(handles=[pc.Patch(color = 'blue',label = 'Pasadena'),pc.Patch(color = 'red',label = 'Palos Verdes')])
    
    #"Main" method
    
    pasCO2Shaved,pasTimeShaved = shave(pasCO2,pasTime)
    pvCO2Shaved,pvTimeShaved = shave(pvCO2,pvTime)
    
    pasCO2removed = doMath(pasCO2Shaved)
    pvCO2removed = doMath(pvCO2Shaved)
    
    pas2004CO2,pas2004Time = removeBefore2004(pasCO2removed,pasTimeShaved)
    pasFirstTime,pasSecondTime = split(pas2004Time)
    pasFirstHalf,pasSecondHalf = split(pas2004CO2)
    pvFirstHalf,pvSecondHalf = split(pvCO2removed)
    
    #seasonal cycles
    pasSeason,pasSemi,pasTrend,pasAll,pasTimeProper = seasonalCycle(pas2004CO2,pas2004Time)
    pvSeason,pvSemi,pvTrend,pvAll,pvTimeProper = seasonalCycle(pvCO2removed,pvTime)
    """
    #logNormal
    plt.figure()
    logNormal(pasFirstHalf,-3,3,1000,"Pasadena First Half",'blue','log','logspace')
    logNormal(pasSecondHalf,-3,3,1000,"Pasadena Second Half",'red','log','logspace')
    logNormal(pvFirstHalf,-3,3,1000,"PV First Half",'green','log','logspace')
    logNormal(pvSecondHalf,-3,3,1000,"PV Second Half",'orange','log','logspace')
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
    """
    """
    #season
    plt.figure()
    plt.plot(pasTimeProper,pasSeason,'blue')
    plt.plot(pvTimeProper,pvSeason,'red')
    seasonalLegend()
    plt.title("PV/Pas Season")
    genericAxes()
    plt.show()
    
    #semi
    plt.figure()
    plt.plot(pasTimeProper,pasSemi,'blue')
    plt.plot(pvTimeProper,pvSemi,'red')
    seasonalLegend()
    plt.title("PV/Pas Semi")
    genericAxes()
    plt.show()
    
    #trend
    plt.figure()
    plt.plot(pasTimeProper,pasTrend,'blue')
    plt.plot(pvTimeProper,pvTrend,'red')
    seasonalLegend()
    plt.title("PV/Pas Trend")
    genericAxes()
    plt.show()
"""
    """
    #all
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
    """
    dFrame = averageYears(pasAll,pasTimeProper)
    dFrame2 = makeConsecutiveDates(dFrame)
    

#-------------------------------------------------La Jolla
elif site == "La Jolla":
    if local == True:
        data_book = pandas.ExcelFile("C:\Users\Kyle W\Desktop\dailyflasks_co2_ljo 07302017.xlsx")
    else:
        data_book = ""
    
    laJolla_data = data_book.parse("data")
    laJolla_time = laJolla_data["Excel date"]
    laJolla_CO2 = laJolla_data["good CO2 (ppm)"]
    
    laJollaCO2removed = doMath(laJolla_CO2)
    
else:
    print "Dead end!"