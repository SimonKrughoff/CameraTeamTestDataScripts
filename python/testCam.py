#!/usr/bin/env python
import numpy, argparse
from lsst.afw.cameraGeom.fitsUtils import DetectorBuilder, setByKey, getByKey
from lsst.ip.isr import AssembleCcdTask
import lsst.afw.math as afwMath
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

class TestCamDetectorBuilder(DetectorBuilder):
    def _sanitizeHeaderMetadata(self, metadata, clobber):
        naxis1 = 544
        naxis2 = 2048
        #all amps have the same biassec
        setByKey('BIASSEC', '[523:544, 1:2002]', metadata, clobber)
        #Get channel number and convert to zero index
        channel = getByKey('CHANNEL', metadata)-1
        if channel is None:
            raise ValueError("Channel keyword not found in header")
        ny = channel//8
        sign = 1 if ny%2 else -1
        nx = 7*ny - sign*(channel%8)
        print channel, nx, ny
        setByKey('DTV1', nx*naxis1, metadata, clobber)
        setByKey('DTV2', ny*naxis2, metadata, clobber)
        #map to the keyword expected for this value
        setByKey('DTM1_1', metadata.get('LTM1_1'), metadata, clobber)
        setByKey('DTM2_2', metadata.get('LTM2_2'), metadata, clobber)
        self._defaultSanitization(metadata, clobber)

class imageSource(object):
    def __init__(self, exposure):
        self.exposure = exposure
        self.image = self.exposure.getMaskedImage().getImage()
    def getCcdImage(self, det, imageFactory, binSize):
        return afwMath.binImage(self.image, binSize)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example to assemble amps from disk")
    parser.add_argument("detectorFile", type=str,
                        help="FITS file to read for detector level metadata")
    parser.add_argument("--detectorFileExt", type=int, default=0,
                        help="Extension to use in reading detector metadata")
    parser.add_argument("--ampFile", type=str, default='',
                        help="FITS file to use for amp metadata and images.  Default is to use the"+\
                        " detectorFile")
    parser.add_argument("--numAmps", type=int, default=16,
                        help="Number of amp segments.  Default is 16")
    parser.add_argument("--doGainCorrection", action='store_true',
                        help="Make a rudimentary per amp gain correction? Default is False.")
    parser.add_argument("--displayTrimmed", action='store_true',
                        help="Display trimmed detector in ds9?")
    parser.add_argument("--displayUnTrimmed", action='store_true',
                        help="Display un-trimmed detector in ds9?")
    args = parser.parse_args()
    imDict = {}
    afilelist = []
    dfilename = '%s[%i]'%(args.detectorFile, args.detectorFileExt)
    if not args.ampFile:
        args.ampFile = args.detectorFile
    for i in range(args.numAmps):
        filename = '%s[%i]'%(args.ampFile, i+1)
        md = afwImage.readMetadata(filename)
        afilelist.append(filename)
        imDict[md.get('EXTNAME')] = afwImage.ImageF(filename)

    db = TestCamDetectorBuilder(dfilename, afilelist, inAmpCoords=True, clobberMetadata=True)
    det = db.buildDetector()
    assembleInput = {}
    for amp in det:
        im = imDict[amp.getName()]
        oscanim = im.Factory(im, amp.getRawHorizontalOverscanBBox())
        oscan = numpy.median(oscanim.getArray())
        imArr = im.getArray()
        imArr -= oscan
        #Calculate and correct for gain
        if args.doGainCorrection:
            # Buffer so edge rolloff doesn't interfere
            medCounts = numpy.median(imArr[30:-30][30:-30])
            stdCounts = numpy.std(imArr[30:-30][30:-30])
            gain = medCounts/stdCounts**2
            imArr *= gain
        assembleInput[amp.getName()] = db.makeExposure(im)

    assembleConfig = AssembleCcdTask.ConfigClass()

    if args.displayTrimmed:
        assembleConfig.doTrim = True
        assembler = AssembleCcdTask(config=assembleConfig)
        resultExp = assembler.assembleCcd(assembleInput)
        camGeomUtils.showCcd(resultExp.getDetector(), imageSource(resultExp), frame=0)

    if args.displayUnTrimmed:
        assembleConfig.doTrim = False
        assembler = AssembleCcdTask(config=assembleConfig)
        resultExp = assembler.assembleCcd(assembleInput)
        camGeomUtils.showCcd(resultExp.getDetector(), imageSource(resultExp), frame=1)


