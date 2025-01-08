#@ String (label="C1", choices={"Nuclei", "PGC Nuclei", "PGC Cytoplasm", "Measure", "None"}, value="DAPI", style="listBox") c1op
#@ String (label="C2", choices={"Nuclei", "PGC Nuclei", "PGC Cytoplasm", "Measure", "None"}, value="PGC Marker", style="listBox") c2op
#@ String (label="C3", choices={"Nuclei", "PGC Nuclei", "PGC Cytoplasm", "Measure", "None"}, value="Measure", style="listBox") c3op
#@ String (label="C4", choices={"Nuclei", "PGC Nuclei", "PGC Cytoplasm", "Measure", "None"}, value="Measure", style="listBox") c4op
#@ Boolean (label="Plot", value=True) showplots
#@ Boolean (label="Histogram Table", value=True) histtable
#@ Boolean (label="Raw Table", value=False) rawtable


import math as maths
from collections import OrderedDict

from ij import IJ, ImagePlus
from ij.plugin.filter import ThresholdToSelection
from ij.process import ImageProcessor, ShortProcessor, FloatProcessor, ColorProcessor, AutoThresholder, Blitter, StackStatistics
from ij.measure import Calibration, ResultsTable
from ij.gui import Roi, ShapeRoi, PolygonRoi, Overlay

from java.awt import Color, Polygon, Font


nucR = 5.0	# estimated nucleus radius (microns)
eR = 2.0 # erode radius (microns)

def plot(histogram, stats, colour):
	nbins = len(histogram[0])
	step = maxv / float(nbins)
	
	width = 512
	halfW = width//2
	height = 512
	pad = 4

	barH = height/nbins

	maxN = max(histogram[1])
	totalN = sum(histogram[1])
	mean = stats.mean

	q1 = None
	median = None
	q3 = None
	cumn = 0
	for i in range(nbins):
		cumn += histogram[1][i]
		if q1 is None and cumn >= totalN*0.25:
			q1 = histogram[0][i]
		if median is None and cumn >= totalN*0.5:
			median = histogram[0][i]
		if q3 is None and cumn >= totalN*0.75:
			q3 = histogram[0][i]
			break
	
	meanY = height - int(mean/maxv * height)
	q1Y = height - int(q1/maxv * height)
	medianY = height - int(median/maxv * height)
	q3Y = height - int(q3/maxv * height)
	
	cp = ColorProcessor(width+2*pad, height+2*pad)
	cp.setColor(Color(192,192,192))
	cp.fillRect(0,0, cp.getWidth(), cp.getHeight())

	left = []
	right = []
	for i in range(nbins):
		y = height - int(i/float(nbins-1) * (height-1)) + pad
		w = int(maths.ceil( (histogram[1][i]/float(maxN)) * (halfW-pad))) # ceil to give 1 px width for small counts, 0 for 0 counts
		left.append([halfW-w+pad,y])
		right.append([halfW+w+pad,y])
	
	# filled area
	poly = Polygon()
	for p in left:
		if p[0] >= halfW+pad-1: continue
		poly.addPoint(p[0], p[1]+barH)
		poly.addPoint(p[0], p[1])
	for p in reversed(right):
		if p[0] <= halfW+pad+1: continue
		poly.addPoint(p[0], p[1])
		poly.addPoint(p[0], p[1]+barH)
	pr = PolygonRoi(poly, PolygonRoi.POLYGON)
	
	cp.setColor(colour)
	cp.fill(pr)
	
	# bars
	'''cp.setColor(colour)
	for i in range(len(left)):
		cp.fillRect(left[i][0],left[i][1], right[i][0]-left[i][0],barH)'''
	
	# IQR box
	cp.setColor(colour.darker())
	cp.fillRect(halfW-64+pad,q3Y+pad, 128,q1Y-q3Y+pad)
	
	# median line
	cp.setColor(Color.BLACK)
	cp.fillRect(halfW-64+pad,medianY+pad, 128,barH)
	
	# mean
	cp.setColor(Color.BLACK)
	cp.fillOval(halfW-barH+pad,meanY-barH+pad, 2*barH+1,2*barH+1)
	
	# CI line
	Z = 1.96 # 0.95 confidence
	ci = Z * (stats.stdDev/maths.sqrt(totalN))
	ciY = int(ci/maxv * height)
	cp.setColor(Color.WHITE)
	cp.drawLine(halfW+pad,meanY+pad-ciY, halfW+pad,meanY+pad+ciY)
	cp.drawLine(halfW+pad-16,meanY+pad-ciY, halfW+pad+16,meanY+pad-ciY)
	cp.drawLine(halfW+pad-16,meanY+pad+ciY, halfW+pad+16,meanY+pad+ciY)
	
	return cp


def getRoi(ip, sigmaPx, method, er):
	proc = ip.duplicate()
	proc.blurGaussian(sigmaPx)
	stats = proc.getStatistics()
	hist = [int(h) for h in stats.getHistogram()]
	thresh = AutoThresholder().getThreshold(method, hist)
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min

	proc.threshold( int(thresh) )
	mask = proc.convertToByte(False)

	for e in range(er):
		mask.dilate() # inverse erode
	for e in range(er):
		mask.erode() # inverse dilate

	mask.setThreshold( 255, 255, ImageProcessor.NO_LUT_UPDATE )
	tts = ThresholdToSelection()
	roi = tts.convert(mask)
	
	return roi


def getValues(ip, roi):
	values = []
	for p in roi.getContainedPoints():
		values.append(ip.getf(p.x,p.y))
	return values


def histogram(values):
	nbins = 256
	step = maxv // (nbins-1)
	binValues = [i*step for i in range(nbins)]
	binCounts = [0 for i in range(nbins)]
	for v in values:
		try:
			bini = int(v/maxv * (nbins-1))
			binCounts[bini] += 1
		except:
			print(v,maxv,bini,nbins)
			exit("nnnnnn")

	totalN = float(sum(binCounts))
	binFreq = [n/totalN for n in binCounts]
	return [binValues, binCounts, binFreq]


ops = [None,c1op,c2op,c3op,c4op]
if ops.count("Nuclei") != 1 or ops.count("PGC Nuclei") != 1 or ops.count("PGC Cytoplasm") > 1:
	IJ.error("One Nuclei and one PGC Nuclei channel should be set, one PGC Cytoplasm channel is optional.")
	exit("invalid channel operations")

nucleiC = ops.index("Nuclei")
pgcC = ops.index("PGC Nuclei")
pgccytoC = ops.index("PGC Cytoplasm")  if "PGC Cytoplasm" in ops else -1
measureCs = [c for c in range(len(ops)) if ops[c]=="Measure"]

imp = IJ.getImage()
stack = imp.getStack()
ol = Overlay()
rt = ResultsTable.getResultsTable()
ht = ResultsTable() if histtable else None
rawrt = ResultsTable() if rawtable else None

W = imp.getWidth()
H = imp.getHeight()
C = imp.getNChannels()
Z = imp.getNSlices()
cal = imp.getCalibration()
unitA = " ("+cal.getUnit()+u"\u00b2"+")"
nucRpx = nucR/cal.pixelWidth
nucApx = maths.pi*nucRpx*nucRpx

ss = StackStatistics(imp)
maxv = 2**12-1
if ss.max>maxv:
	maxv = 2**16-1

images = {}
for c in range(1,C+1):
	if Z == 1:
		images[c] = stack.getProcessor(c)
	else:
		proj = ShortProcessor(W,H)
		for z in range(1,Z):
			proj.copyBits(stack.getProcessor(imp.getStackIndex(c,z,1)), 0,0, Blitter.MAX)
		images[c] = proj

nucleiRoi = getRoi(images[nucleiC], 1, AutoThresholder.Method.Triangle, 0)
if nucleiRoi is None:
	exit("no nuclei")

epx = int(eR/cal.pixelWidth)
pgcRoi = getRoi(images[pgcC], 1, AutoThresholder.Method.MaxEntropy, epx)
pgccytoRoi = getRoi(images[pgccytoC], 2, AutoThresholder.Method.Triangle, 0) if pgccytoC > -1 else None

if pgccytoRoi != None:
	pgcRoi = ShapeRoi(pgcRoi).and(ShapeRoi(pgccytoRoi))

nonpgcRoi = ShapeRoi(nucleiRoi).not(ShapeRoi(pgcRoi))

nonpgcRoi.setStrokeColor(Color.YELLOW)
ol.add(nonpgcRoi)
#nucleiRoi.setStrokeColor(Color.CYAN)
#ol.add(nucleiRoi)
if pgccytoRoi is not None:
	pgccytoRoi.setStrokeColor(Color.CYAN)
	ol.add(pgccytoRoi)
pgcRoi.setStrokeColor(Color.MAGENTA)
ol.add(pgcRoi)

plots = OrderedDict()
row = rt.getCounter()
for c in measureCs:
	images[c].setRoi(nucleiRoi)
	nucleiStats = images[c].getStatistics()
	images[c].setRoi(pgcRoi)
	pgcStats = images[c].getStatistics()
	pgcValues = getValues(images[c], pgcRoi)
	images[c].setRoi(nonpgcRoi)
	nonpgcStats = images[c].getStatistics()
	nonpgcValues = getValues(images[c], nonpgcRoi)
	if c == measureCs[0]:
		rt.setValue("Image", row, imp.getTitle())
		rt.setValue("Channel", row, c)
		rt.setValue("Nuclei Area"+unitA, row, nucleiStats.area*cal.pixelWidth*cal.pixelHeight)
		rt.setValue("PGC Nuclei Area"+unitA, row, pgcStats.area*cal.pixelWidth*cal.pixelHeight)
		rt.setValue("Non-PGC Area"+unitA, row, nonpgcStats.area*cal.pixelWidth*cal.pixelHeight)
		rt.setValue("Nuclei", row, nucleiStats.area/nucApx)
		rt.setValue("PGC Nuclei", row, pgcStats.area/nucApx)
		rt.setValue("Non-PGCs", row, nonpgcStats.area/nucApx)
	rt.setValue("C"+str(c)+" Nuclei Mean", row, nucleiStats.mean)
	rt.setValue("C"+str(c)+" PGC Nuclei Mean", row, pgcStats.mean)
	rt.setValue("C"+str(c)+" Non-PGC Mean", row, nonpgcStats.mean)

	if showplots or histtable:
		pgcHist = histogram(pgcValues)
		nonpgcHist = histogram(nonpgcValues)
	
	if showplots:
		plots["C"+str(c)+" PGC Nuclei"] = plot(pgcHist, pgcStats, Color.MAGENTA)
	 	plots["C"+str(c)+" Non-PGC Nuclei"] = plot(nonpgcHist, nonpgcStats, Color.YELLOW)
	if histtable:
		for i in range(256):
			ht.setValue("Bin",i,i)
			ht.setValue("Intensity",i,pgcHist[0][i])
			ht.setValue("C"+str(c)+" PGC Nuclei n Pixels",i,pgcHist[1][i])
			ht.setValue("C"+str(c)+" PGC Nuclei f",i,pgcHist[2][i])
			ht.setValue("C"+str(c)+" Non-PGC Nuclei n Pixels",i,nonpgcHist[1][i])
			ht.setValue("C"+str(c)+" Non-PGC Nuclei f",i,nonpgcHist[2][i])
	if rawtable:
		maxn = max(len(pgcValues), len(nonpgcValues))
		for i in range(maxn):
			rawrt.setValue("PGC Nuclei Pixel Values", i, pgcValues[i] if i < len(pgcValues) else "")
			rawrt.setValue("Non-PGC Nuclei Pixel Values", i, nonpgcValues[i] if i < len(nonpgcValues) else "")

imp.setOverlay(ol)
rt.show("Results")
if histtable:
	ht.show("PGC Measure "+imp.getTitle())
	
if rawtable:
	rawrt.show("PGC Measure Raw Values "+imp.getTitle())

if showplots:
	head = 32
	pad = 20
	width = plots[plots.keys()[0]].getWidth()
	height = plots[plots.keys()[0]].getHeight()
	pcp = ColorProcessor(sum(plots[key].getWidth()+pad for key in plots.keys())+pad, head+height+pad)
	sumW = pad
	pcp.setFont(Font(Font.SANS_SERIF, Font.PLAIN, 12))
	for key in plots.keys():
		pcp.setColor(Color.WHITE)
		pcp.drawString(key, sumW+pad, pad)
		pcp.drawString("0", sumW+width/2-3, head+height+14)
		pcp.drawString("%.i"%maxv, sumW+width/2-20, head)
		
		plot = plots.get(key)
		pcp.copyBits(plot, sumW,head, Blitter.COPY)
		sumW += plot.getWidth()+pad
	
	ImagePlus("PGC Measure "+imp.getTitle(), pcp).show()
