import os
import sys
import optparse
import ROOT
from array import array

NLAYERS=28

def fixExtremities(h,addOverflow=True,addUnderflow=True):
    """increments the first and the last bin to show the under- and over-flows"""
    if addUnderflow :
        fbin  = h.GetBinContent(0) + h.GetBinContent(1)
        fbine = ROOT.TMath.Sqrt(h.GetBinError(0)*h.GetBinError(0) + h.GetBinError(1)*h.GetBinError(1))
        h.SetBinContent(1,fbin)
        h.SetBinError(1,fbine)
        h.SetBinContent(0,0)
        h.SetBinError(0,0)
    if addOverflow:
        nbins = h.GetNbinsX();
        fbin  = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
        fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)*h.GetBinError(nbins)  + h.GetBinError(nbins+1)*h.GetBinError(nbins+1))
        h.SetBinContent(nbins,fbin)
        h.SetBinError(nbins,fbine)
        h.SetBinContent(nbins+1,0)
        h.SetBinError(nbins+1,0)

def saveProfilesFor(distName,inUrl, thrScan, outDir, quantile=0.5):
    """analyzes the histograms and dumps a small tree for further analysis"""

    #prepare output
    outUrl=os.path.join(outDir,'%s_summary_%s'%(distName,os.path.basename(inUrl)))
    fOut=ROOT.TFile.Open(outUrl,'RECREATE')
    fOut.cd()
    tree = ROOT.TNtuple(distName,distName,"layer:ieta:eta:thr:median:medianUnc")

    print 'Analyzing',inUrl,'output in',outUrl

    #profile the input
    fIn=ROOT.TFile.Open(inUrl)
    for thr in thrScan:
        print '\t threshold',thr

        fdir=fOut.mkdir('thr%d'%thr)
        for ilay in xrange(1,NLAYERS+1):

            distH=fIn.Get("ana/sd_EE_layer%d_thr%d_%s"%(ilay,thr,distName))

            for ybin in xrange(1,distH.GetNbinsY()+1):

                #project and compute quantiles
                fIn.cd()
                etaMin,etaMax=distH.GetYaxis().GetBinLowEdge(ybin),distH.GetYaxis().GetBinUpEdge(ybin)
                projName="%s_ithr%d_ieta%d_profile"%(distName,thr,ybin)
                px=distH.ProjectionX(projName,ybin,ybin);
                px.SetLineWidth(2);
                px.SetTitle( "Layer %d, %3.2f<|#eta|<%3.2f"%(ilay,etaMin,etaMax) )
                fixExtremities(px)
                xq = array('d', [0.5])
                yq = array('d', [0.])
                px.GetQuantiles(1,yq,xq)

                #save to output
                fOut.cd()
                tree.Fill(ilay,ybin,0.5*(etaMin+etaMax),thr,yq[0],1.253*px.GetMeanError())
                fdir.cd()
                px.SetDirectory(fdir)
                px.Write()

    fIn.Close()
    fOut.cd()
    tree.Write()
    fOut.Close()

def saveBusyProfilesFor(inUrl,outDir):
    """analyzes the histograms and dumps a small tree for further analysis"""

    #prepare output
    outUrl=os.path.join(outDir,'busy_summary_%s'%(os.path.basename(inUrl)))
    fOut=ROOT.TFile.Open(outUrl,'RECREATE')
    fOut.cd()
    tree = ROOT.TNtuple('busy','busy',"layer:ieta:eta:thr:median:medianUnc")

    print 'Analyzing',inUrl,'output in',outUrl

    #profile the input
    fIn=ROOT.TFile.Open(inUrl)
    for ilay in xrange(1,NLAYERS+1):

            totalch=fIn.Get("ana/sd_EE_layer%d_nch"%ilay)
            distH=fIn.Get("ana/sd_EE_layer%d_busy"%ilay)
            for ybin in xrange(1,distH.GetNbinsY()+1):
                fIn.cd()
                etaMin,etaMax=distH.GetYaxis().GetBinLowEdge(ybin),distH.GetYaxis().GetBinUpEdge(ybin)
                px=distH.ProjectionX("px",ybin,ybin)
                busy,busyUnc = px.GetBinContent(2),px.GetBinError(2)
                chcount      = totalch.GetBinContent(ybin)
                fbusy        = busy/chcount                
                fbusyUnc     = busyUnc/chcount
                px.Delete()
                fOut.cd()
                tree.Fill(ilay,ybin,0.5*(etaMin+etaMax),0,fbusy,fbusyUnc)

    fIn.Close()
    fOut.cd()
    tree.Write()
    fOut.Close()


def main():
    """parse inputs and run the analysis"""

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          dest='input',       help='input file [%default]',               default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory [%default]',         default='summary',  type='string')
    parser.add_option(      '--var',         dest='var',         help='variable to profile  [%default]',     default='occ',      type='string')
    parser.add_option(      '--thrScan',     dest='thrScan',     help='thresholds to scan (csv) [%default]', default='2,4,20', type='string')
    parser.add_option(      '--doBusy',      dest='doBusy',      help='do busy profile [%default]',          default=False,      action='store_true')
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s'%opt.output)
    saveProfilesFor(distName=opt.var,inUrl=opt.input,
                    thrScan=[int(x) for x in opt.thrScan.split(',')],
                    outDir=opt.output,
                    quantile=0.5)
    if opt.doBusy:
        saveBusyProfilesFor(inUrl=opt.input,outDir=opt.output)


if __name__ == "__main__":
    sys.exit(main())
