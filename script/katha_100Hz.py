import numpy as np
import struct
from decimal import Decimal
import os
import matplotlib.pyplot as plt
from statistics import mean
from math import sqrt
import shutil
# Set path to folder where your .TXT files are located here
# If script is where your .TXT are located then leave it as filePath = r'.'
filePath = r'.'
#filePath = 'C:/Users/745899/Documents/Islington/BDS_Sensors_13_02_2023'


def main():
    txtFiles = getTXTFilesFromDirectory(filePath)
    createDirectories()
    for i in range(len(txtFiles)): # Repeats all actions until all .TXT files have been read
        print("Processing file {} of {}".format(i+1, len(txtFiles)))
        fileName = txtFiles[i]
        listOfValues = get100HzBinaryFileValues(fileName)
        if listOfValues:
            exportToCSV(listOfValues, fileName)
            exportToPDF(listOfValues, fileName)
            moveTXTtoTXT(fileName)
        else:
            print("     250Hz or Invalid/Unreadable file {}".format(fileName))


# Finds all .TXT files from the input directory 'filePath'
def getTXTFilesFromDirectory(path):
    txtFiles = []
    for file in os.listdir(path):
        if file.endswith('.txt'):
            txtFiles.append(file)
    return txtFiles

# Opens binary .TXT files and gets all their values
def get100HzBinaryFileValues(fileName):
    fmt = '2B H L 22f 4B'
    listOfValues = []
    fileLoc = os.path.join(filePath, fileName)
    with open(fileLoc, 'rb') as f: # Open binary file
        while True:
            row = f.read(2 + 2 + 4 + 22 * 4 + 4) # Number of bytes on one line
            if not row: # Check for EOF
                break
            try:
                # Read values into variables and input to a matrix
                rand1, rand2, sampleRate, timestamp, *float_numbers, calMag, calAcc, calGyro, calImu = struct.unpack(fmt, row)
                if (sampleRate == 100): # Only continues if the binary files are 100Hz
                    listOfValues.extend([[sampleRate, timestamp, *float_numbers, calMag, calAcc, calGyro, calImu]])
                else:
                    break
            except Exception as e:
                error = e
    return listOfValues

# Gets binary file values and inserts into a .CSV file
def exportToCSV(listOfValues, fileName):
    listLength = len(listOfValues)
    for i in range(listLength):
        listOfValues[i].pop(0) # Removes sampleRate from list (not used in .CSV)
        tm = (listOfValues[i][0] - listOfValues[0][0]) / (60 * 1000)
        ts = tm * 60
        if i>=1:
            listOfValues[i][0] = ts
    listOfValues[0][0] = 0
    csvHeader = "Time [s],PL [mbar],TL [C],PC [mbar],TC [C],PR [mbar],TR [C],EX [deg],EY [deg],EZ [deg],QW [-],QX [-],QY [-],QZ [-],MX [microT],MY [microT],MZ [microT],AX [m/s2],AY [m/s2],AZ [m/s2],RX [rad/s],RY [rad/s],RZ [rad/s],CSM,CSA,CSR,CSTOT"
    csvFmt = "%0.2f, %0.1f, %0.2f, %0.1f, %0.2f, %0.1f, %0.1f, %0.2f, %0.3f, %0.2f, %0.5f, %0.7f, %0.5f, %0.5f, %0.3f, %0.3f, %0.3f, %0.2f, %0.2f, %0.2f, %0.7f, %0.7f, %0.7f, %d, %d, %d, %d"
    fileName = fileName.rsplit(".",2)
    fileName[0] += ".csv"
    exportPath = os.path.join(filePath, "CSV100Hz\\")
    exportPath += fileName[0]
    try:
        np.savetxt(exportPath, listOfValues, fmt=csvFmt, delimiter=",", comments='', header=csvHeader)
        print ("Done")
    except Exception as e:
        print(e)

# Creates the folders "CSV100Hz", "TXT100Hz" and "PDF"
def createDirectories():
    path = ["CSV100Hz", "TXT100Hz", "PDF"]
    for i in range(len(path)):
        out = os.path.join(filePath,path[i])
        if not os.path.exists(out):
            os.makedirs(out)

# Moves 100Hz binary .TXT files to the folder "TXT100Hz" after processing
def moveTXTtoTXT(fileName):
    fileLoc = os.path.join(filePath, fileName)
    shutil.move(fileLoc, "{}\TXT100Hz\{}".format(filePath,fileName))

# Creates two plots and exports them to a A4 sized .PDF file
def exportToPDF(listOfValues, fileName):
    listLength = len(listOfValues)
    ts = [] # Timestamp
    pl = [] # PL [mbar]
    pc = [] # PC [mbar]
    pr = [] # PR [mbar]
    mn = [] # mean
    for i in range(listLength):
        ts.append(listOfValues[i][0])
        pl.append(listOfValues[i][1])
        pc.append(listOfValues[i][3])
        pr.append(listOfValues[i][5])
        mn.append(mean([pc[i],pc[i],pr[i]]))
    fig = plt.figure(figsize=(8.3,11.69))
    plt.rcParams.update({"font.size": 10})
    plot1 = fig.add_subplot(2,1,1)
    # "rasterized=True" makes smaller PDFs, otherwise takes ~5 seconds to load the PDF
    plot1.plot(ts, pl, "-r", label="Left", rasterized=True) 
    plot1.plot(ts, pc, "-b", label="Center", rasterized=True)
    plot1.plot(ts, pr, "-k", label="Right", rasterized=True)
    plot1.plot(ts, mn, ".", color="#00FE00", markersize = "2", label="Mean", rasterized=True)
    
    plot1.set_xlabel("Time (s)")
    plot1.set_ylabel("Total Pressure (mbar)")
    plot1.grid()
    plot1.set_title(fileName, fontsize=10, fontweight="bold")
    plot1.legend()

    tick = 0
    xticks = []
    while tick < ts[len(ts)-1]: # Setting length of x and y axis
        xticks.append(tick)
        tick = tick + 100
    else:
        xticks.append(ts[len(ts)-1])
        plot1.set_xticks(xticks) 
        plot1.set_yticks(range(500,int(max(pl)+500),500))

    plot2 = fig.add_subplot(2,1,2)
    acc_mag = []
    for i in range(len(ts)):
        acc_mag.append(sqrt(pow(listOfValues[i][17],2) + pow(listOfValues[i][18],2) + pow(listOfValues[i][19],2)))
    plot2.plot(ts, acc_mag, "-r", label="Acc XYZ", rasterized=True)
    plot2.set_xlabel("Time (s)")
    plot2.set_ylabel("Accel Mag (m/s\u00b2)")
    plot2.grid()
    plot2.legend()

    plot2.set_xticks(xticks) 
    plot2.set_yticks(range(0,int(max(acc_mag))+10,10))

    fileName = fileName.rsplit(".",2)
    fileName[0] += ".pdf"
    exportPath = os.path.join(filePath, "PDF\\")
    exportPath += fileName[0]
    fig.subplots_adjust(hspace=0.5)
    plt.savefig(exportPath, dpi=300) # Export as .PDF file
    #plt.show() # Uncomment to show graphs

if __name__ == "__main__":
    main()