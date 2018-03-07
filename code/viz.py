# Visualization tools for graphrat simulator
# For printing ASCII representation of maze 
# And for showing as a heat map
import math
import sys
import curses
import time
import datetime

# Some installations don't support Tkinter and PIL libraries.
# Import them only if needed

specialImported = False

def importSpecial():
    global specialImported, Tkinter, Image, ImageDraw
    if not specialImported:
        import Tkinter
        from PIL import Image, ImageDraw
    specialImported = True

class VizMode:
    nothing, ascii, heatmap, both, error = range(5)

    def parse(self, name):
        if len(name) != 1:
            return self.error
        elif name == 'n':
            return self.nothing
        elif name == 'a':
            return self.ascii
        elif name == 'h':
            return self.heatmap
        elif name == 'b':
            return self.both
        else:
            return self.error

    def doHeatMap(self, v):
        return v in [self.heatmap, self.both]

    def doASCII(self, v):
        return v in [self.ascii, self.both]

class Formatter:

    k = 0
    digits = 1
    separator = ""
    file = None
    x = 0
    y = 0
    stdscr = None
    tlast = None
    display = None

    def __init__(self, k, maxval, fname = "", viz = VizMode.both):
        self.k = k
        self.digits = int(math.ceil(math.log10(maxval+1)))
        # Pattern that separates lines
        s = "-" * self.digits + "+"
        self.separator = "+" + s * self.k
        self.file = None
        self.stdscr = None
        self.display = None
        vm = VizMode()
        if fname == "":
            self.file = sys.stdout 
        else:
            try:
                self.file = open(fname, "w")
            except:
                self.file = None
                print "Could not open file '%s'" % fname
                return
        if vm.doASCII(viz):
            self.stdscr = curses.initscr()
            curses.noecho()
            curses.cbreak()
        self.tlast = datetime.datetime.now()
        if vm.doHeatMap(viz):
            self.display = Display(k = self.k, maxval = maxval)
        self.reset()

    def reset(self):
        self.x = 0
        self.y = 0

    def finish(self, fname = ""):
        if self.file is not None and self.file not in [sys.stdout, sys.stderr]:
            self.file.close()
            self.file = None
        if fname != "" and self.display is not None:
            self.display.capture(fname)
        self.finishDynamic(wait = False)

    # Complete dynamic part of visualization
    def finishDynamic(self, wait = False):
        if self.stdscr is not None:
            if wait:
                self.printLine("[Hit any key to exit]")
                self.stdscr.getch()
            curses.nocbreak()
            curses.echo()
            curses.endwin()
        self.stdscr = None

    def printLine(self, s, addReturn = True):
        ps = s
        ss = s
        if s[-1] == '\n':
            ss = s[:-1]
            addReturn = True
        else:
            ps = ps + '\n' if addReturn else ps
        if self.stdscr is not None:
            self.stdscr.addstr(self.y, self.x, ss)
        elif (self.file is not None) and (self.display is None):
            self.file.write(ps)
        if addReturn:
            self.x = 0
            self.y += 1

    def show(self, populations, period = 0.0):
        self.printLine(self.separator)
        for row in range(self.k):
            rstart = row * self.k
            rvals = populations[rstart:rstart+self.k]
            line = "|"
            for v in rvals:
                sdig = "" if v == 0 else"%d" % v
                slen = self.digits - len(sdig)
                llen = (slen+1)/2
                rlen = slen/2
                line += " " * llen
                line += sdig
                line += " " * rlen
                line += "|"
            self.printLine(line)
            self.printLine(self.separator)
        if self.stdscr is not None:
            self.stdscr.refresh()
        if period > 0:
            if self.tlast is not None:
                tnow = datetime.datetime.now()
                d = tnow - self.tlast
                dsecs = 24 * 3600 * d.days + d.seconds + 1e-6 * d.microseconds
                secs = period - dsecs
                if secs > 0:
                    time.sleep(secs)
        if self.display:
            self.display.setColors(populations)
        self.tlast = datetime.datetime.now()

# For generating heatmap

def cstring(c):
    fields = [("%.2x" % int(255*x)) for x in c]
    return "#" + "".join(fields)

class Colors:
    black = (0.0, 0.0, 0.0)
    blue =  (0.0, 0.0, 1.0)
    green = (0.0, 1.0, 0.0)
    cyan =  (0.0, 1.0, 1.0)
    red =   (1.0, 0.0, 0.0)
    magenta = (1.0, 0.0, 1.0)
    yellow = (1.0, 1.0, 0.0)
    white = (1.0, 1.0, 1.0)
    lightred = (1.0, 0.5, 0.5)


class HeatMap:
    maxval = 1
    logmax = 0
    colorList = (Colors.magenta, Colors.blue, Colors.cyan, Colors.green, Colors.yellow, Colors.red)
    weightList = (0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    # Below splitVal, interpolate over color range.
    # Above, interpolate between final two colors
    splitVal = 0.5

    def __init__(self, maxval = 100, avgval = 10):
        self.scaleList = []
        for i in range(len(self.colorList)):
            self.scaleList.append(map(lambda (c): c * self.weightList[i], self.colorList[i]))
        self.maxval = maxval
        self.logmax = math.log(self.maxval+1)
        self.splitval = math.log(2.0 * avgval)/self.logmax
        
    # Perform interpolation of number in range [0, 1.0) to get color
    def interpolate(self, x):
        x = min(x, 1.0)
        x = max(x, 0.0)
        upper = x >= self.splitVal
        if upper:
            lcolor = self.scaleList[-2]
            rcolor = self.scaleList[-1]
            point = (x - self.splitVal) / (1.0 - self.splitVal)
        else:
            x = x / self.splitVal
            segs = len(self.scaleList) - 2
            sx = x * segs
            interval = int(sx)
            point = sx - interval
            lcolor = self.scaleList[interval]
            rcolor = self.scaleList[interval+1]
        color = [rcolor[idx] * point + lcolor[idx] * (1-point) for idx in range(3)]
        return cstring(color)

    def genColor(self, val):
        if val == 0:
            return cstring(Colors.black)
        scale = math.log(val) / self.logmax
        return self.interpolate(scale)

    def genColors(self, vlist):
        clist = [self.genColor(v) for v in vlist]
        return clist

class Display:
    k = 10
    squareSize = 8
    display = None  # TK Window
    frame = None    # Frame within window
    canvas = None   # Canvas within frame
    squares = []    # Set of rectangles, k*k total
    colorList = []  # Most recent set of colors
    hmap = None

    def __init__(self, k = 10,  maxdim = 800, maxval = 100):
        importSpecial()
        self.k = k
        nodeCount = k * k
        loadFactor = maxval / nodeCount
        self.squareSize = maxdim / self.k
        self.display = Tkinter.Tk()
        self.display.title('GraphRat Simulation of %d X %d maze (load factor = %d)' % (k, k, loadFactor))
        self.frame = Tkinter.Frame(self.display)
        self.frame.pack(fill=Tkinter.BOTH)
        self.canvas = Tkinter.Canvas(self.frame, width = k * self.squareSize, height = k * self.squareSize)
        self.canvas.pack(fill=Tkinter.BOTH)
        self.squares = []
        for r in range(0, self.k):
            for c in range(0, self.k):
                (x, y) = self.xyPos(r, c)
                sq = self.canvas.create_rectangle(x, y, x+self.squareSize, y+self.squareSize, width = 0,
                                                  fill = cstring(Colors.black))
                self.squares.append(sq)
        self.hmap = HeatMap(maxval = maxval, avgval = maxval/nodeCount)
        self.update()

    def update(self):
        self.canvas.update()

    def xyPos(self, r, c):
        x = self.squareSize * c
        y = self.squareSize * r
        return (x, y)

    def rowCol(self, idx):
        r = idx / self.k
        c = idx % self.k
        return (r, c)

    def colorSquare(self, idx, color):
        if idx >= 0 and idx < len(self.squares):
            square =  self.squares[idx]
            self.canvas.itemconfig(square, fill = color)

    # Set colors based on counts for each square
    def setColors(self, vlist = []):
        clist = self.hmap.genColors(vlist)
        self.colorList = clist
        for idx in range(len(vlist)):
            self.colorSquare(idx, clist[idx])
        self.update()
            
    def capture(self, fname):
        img = Image.new('RGB', (self.k * self.squareSize, self.k * self.squareSize), "black")
        dimg = ImageDraw.Draw(img)
        for idx in range(len(self.colorList)):
            r, c = self.rowCol(idx)
            x, y = self.xyPos(r, c)
            dimg.rectangle((x, y, x + self.squareSize, y + self.squareSize), fill = self.colorList[idx])
            try:
                img.save(fname)
            except Exception as e:
                print "Could not save image to file %s.  %s" % (fname, e)

    def finish(self):
        self.display.destroy()



